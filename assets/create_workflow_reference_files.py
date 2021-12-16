#! venv/bin/python
# import statements, alphabetic order of main package.
import argparse
from configparser import ConfigParser
from difflib import unified_diff
from errno import ENOENT as errno_ENOENT
from io import IOBase
import json
from os import strerror as os_strerror
import pathlib
import requests
import sys
from warnings import warn as warnings_warn

# third party libraries alphabetic order of main package.
from deepdiff import DeepDiff
from natsort import index_natsorted
from numpy import argsort as np_argsort
import pandas as pd
import yaml

def non_empty_existing_file(file):
    path_file = pathlib.Path(file)
    print(file)
    print(path_file.stat().st_size)
    print("ok")
    print(path_file.stat())
    if not path_file.is_file() and not path_file.is_dir():
        raise FileNotFoundError(errno_ENOENT, os_strerror(errno_ENOENT), file)
    elif not path_file.is_dir() and not path_file.stat().st_size:
        raise OSError("File is empty.")
    elif path_file.is_dir() and not path_file.is_absolute():
        raise OSError("Filepath is expected to be absolute.")
    if path_file.is_dir():
        return file
    else:
        return open(file)

def parse_arguments_and_check(args_in):
    parser = argparse.ArgumentParser(
        description="Migrate translation table to a bed file with sites of interest " +
        "and a yaml file with genotype - phenotype translation."
    )
    parser.add_argument(
        "-b", "--bed", type=non_empty_existing_file, required=False,
        help="Previous bed file. If provided, differences between previous and generated bed file are shown."
    )
    parser.add_argument(
        "-c", "--config_file", 
        type=non_empty_existing_file, 
        default=str(pathlib.Path(__file__).parent) + "/create_workflow_reference_files.ini",
        help="Filepath to INI-config file."
    )
    parser.add_argument(
        "-s", "--config_section", type=str, default="DEFAULT",
        help="The INI-config section used."
    )
    parser.add_argument(
        "-o", "--output_path", type=non_empty_existing_file, required=False, default=pathlib.Path(__file__).cwd(),
        help="Filepath where output is placed. (default: %(default)s)"
    )
    parser.add_argument(
        "-p", "--output_prefix", type=str, required=False,
        help="Output prefix to use as filename in output. Default is filename prefix of translation table."
    )
    parser.add_argument(
        "-y", "--yaml", type=non_empty_existing_file, required=False,
        help="Previous yaml file. If provided, differences between previous and generated yaml file are shown."
    )
    parser.add_argument(
        "translation_table", type=non_empty_existing_file,
        help="Filepath to translaion table file."
    )
    args = parser.parse_args(args_in)
    return(args)


def read_config_section_and_check(section, config_file):
    config_parser = ConfigParser(converters={"jsonloads": json.loads})
    config_parser.read_file(config_file)
    config_section = config_parser[section]
    required_keys = ["ensembl_url", "species"]
    for req_key in required_keys:
        if req_key not in config_section:
            raise KeyError("Required key {} not in config file.".format(req_key))
    return(config_section)


def morph_translation_file(csv_file, dict_rename_cols=None):
    df_translation = pd.read_csv(csv_file)
    if dict_rename_cols:
        df_translation = df_translation.rename(columns=dict_rename_cols)
    required_cols = {"genotype_id", "gene_and_rs_id", "variant_genotype", "phenotype_id", "phenotype_name"}
    if required_cols - set(df_translation.columns):
        raise ValueError("Required columns are missing in translation file: {}".format(
                required_cols - set(df_translation.columns)
            )
        )
    # remove sections
    if any(df_translation[df_translation['genotype_id'].astype(str).str.match('---')]):
        df_translation = df_translation[~df_translation['genotype_id'].astype(str).str.match('---')]
    if any(df_translation.variant_genotype.isna()):
        raise ValueError("Variant_genotype value is required in translation file. Offending rows:\n{}".format(
                df_translation[df_translation.variant_genotype.isna()]
            )
        )
    if any(df_translation.gene_and_rs_id.str.count("_") != 1):
        raise ValueError("Expected separator _ in translation file, column 'gene_and_rs_id'. Offending rows:\n{}".format(
                df_translation[df_translation.gene_and_rs_id.str.count("_") != 1]
            )
        )
    df_translation[['gene', 'rs_id']] = df_translation.gene_and_rs_id.str.split("_", expand=True,)
    if any(
        (df_translation.variant_genotype.str.count(":|/") != 1)
        & (~df_translation.variant_genotype.str.match("Missing", na=False))
        & (df_translation.rs_id.str.startswith("rs", na=False))
    ):
        raise ValueError(
            "Expected single separator ':' or '/' in translation file, column 'variant_genotype'. Offending rows:\n{}".format(
                df_translation[df_translation.variant_genotype.str.count(":|/") != 1]
            )
        )
    if any(
        (df_translation.variant_genotype.str.replace("A|T|C|G|:|-|\\.|\\/", "", regex=True).str.len() != 0)
        & (~df_translation.variant_genotype.str.match("Missing", na=False))
        & (df_translation.rs_id.str.startswith("rs", na=False))
    ):
        raise ValueError(
            "Invalid character in variant genotype. Supported: 'A', 'T', 'C', 'G', ':', '-', '.', '/'"
        )

    df_translation['variant_genotype'] = df_translation.variant_genotype.str.replace(':', "/")
    # split translation table into phenotype and genotype metadata.
    phenotype_cols = [
        "genotype_id", "gene_genotype", "gene", "genotype_realname", "phenotype_id", "phenotype_name"
    ]
    df_phenotypes = (
        df_translation[df_translation.columns.intersection(phenotype_cols)]
        .dropna(axis=0, subset=["gene", "phenotype_name"])
        .set_index("genotype_id", drop=False)
    )
    removed_pt = set(df_translation.genotype_id.tolist()) - set(df_phenotypes.index.tolist())
    if removed_pt:
        print("Following genotype IDs (n={len}) are removed, missing phenotypes or genes. {ids}".format(
            len=len(removed_pt), ids=removed_pt))
    df_genotypes = df_translation.loc[
        df_translation.genotype_id.isin(df_phenotypes.genotype_id.to_list()),
        ["genotype_id", "gene_and_rs_id", "gene", 'rs_id', 'variant_genotype']
    ]
    return(df_phenotypes, df_genotypes)


def get_ensembl_request_response(server, ext, json=None, method="get"):
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    if method == "get":
        r = requests.get(server+ext, headers=headers, json=json)
    elif method == "post":
        r = requests.post(server+ext, headers=headers, json=json)
    if not r.ok:
        r.raise_for_status()
        sys.exit(1)
    return(r.json())


def get_variant_metadata_ensembl(server, ids, species):
    ens_variants = get_ensembl_request_response(
        server=server,
        ext="/variation/{species}".format(species=species),
        json={"ids": ids},
        method="post",
    )
    if not ens_variants:
        warnings_warn("Variation identifiers not found in Ensembl:\n{}".format(ids))
    return(ens_variants)


def get_gene_metadata_ensembl(server, ens_rs_id, location, species):
    client_out = get_ensembl_request_response(
        server=server,
        ext="/overlap/region/{species}/{loc}?feature=gene;logic_name=ensembl_havana_gene_homo_sapiens".format(
            species=species,
            loc=location
        )
    )
    if not client_out:
        warnings_warn("Gene not found for rs ID {} in Ensembl".format(ens_rs_id))
        return(None, None, None)

    gene_id = client_out[0].get("gene_id", None)
    gene_name = client_out[0].get("external_name", None)
    strand = client_out[0].get("strand", None)
    return(gene_name, gene_id, strand)


def get_linked_gene_and_id(df_translation, ens_rs_id, retrieved_rs_synonyms):
    if not df_translation.rs_id.isin([ens_rs_id]).any():
        rs_id_intxn = list(set(df_translation.rs_id) & set(retrieved_rs_synonyms))
        linked_gene = ";".join(df_translation.loc[df_translation.rs_id.isin(rs_id_intxn)].gene.tolist())
        gene_and_rs_id = ";".join(df_translation.loc[df_translation.rs_id.isin(rs_id_intxn)].gene_and_rs_id.tolist())
        warnings_warn(
            "Synonym ID {syn} is used instead of the original rs ID {rs_id} ".format(syn=rs_id_intxn, rs_id=ens_rs_id)
        )
    else:
        linked_gene = df_translation.loc[df_translation.rs_id == ens_rs_id].gene.item()
        gene_and_rs_id = df_translation.loc[df_translation.rs_id == ens_rs_id].gene_and_rs_id.item()
    return(linked_gene, gene_and_rs_id)


def get_data_ensembl(df_translation, ensembl_url, species):
    ens_variants = get_variant_metadata_ensembl(server=ensembl_url, ids=df_translation.rs_id.tolist(), species=species)
    cols = ['chrom', 'start', 'end', 'name', 'gene', 'retrieved_gene', 'retrieved_id', 'ref', 'alt']
    df_ens_metadata = pd.DataFrame(columns=cols)
    for ens_rs_id in ens_variants.keys():
        map_info = ens_variants[ens_rs_id]["mappings"][0]
        gene_name, gene_id, strand = get_gene_metadata_ensembl(
            server=ensembl_url,
            ens_rs_id=ens_rs_id,
            location=map_info["location"],
            species=species
        )
        linked_gene, gene_and_rs_id = get_linked_gene_and_id(
            df_translation=df_translation,
            ens_rs_id=ens_rs_id,
            retrieved_rs_synonyms=ens_variants[ens_rs_id]["synonyms"]
            )
        df_ens_metadata = df_ens_metadata.append(
            {
                'chrom': "chr" + map_info["seq_region_name"],
                'start': map_info["start"]-1,  # transform to zero based.
                'end': map_info["end"],
                'name': gene_and_rs_id,
                'rs_id': ens_rs_id,
                'gene': linked_gene,
                'retrieved_gene': gene_name,
                'retrieved_id': gene_id,
                'ref': map_info["allele_string"].split("/")[0],
                'alt': map_info["allele_string"].split("/", 1)[1],
                'strand': strand,
            }, ignore_index=True
        )
    df_metadata_sort = df_ens_metadata.sort_values(by="chrom", key=lambda x: np_argsort(
        index_natsorted(zip(df_ens_metadata.chrom, df_ens_metadata.start))))
    return(df_metadata_sort)


def get_invalid_genes_and_warn(df_ens_metadata, genes_regex=None):
    lst_filter_gene_names = []
    if genes_regex:
        lst_filter_gene_names += df_ens_metadata.loc[df_ens_metadata.gene.str.contains(genes_regex), "gene"].to_list()
    df_grouped = df_ens_metadata.groupby("gene")['chrom'].nunique()
    if any(df_grouped > 1):
        lst_filter_gene_names += df_grouped.where(df_grouped > 1).dropna().keys().tolist()
        warnings_warn("At least one gene is linked to variants from different chromosomes.\n{counts}".format(
                counts=df_ens_metadata.groupby("gene", as_index=False)["chrom"].nunique()
            )
        )
    df_genes_not_match = df_ens_metadata.query("gene != retrieved_gene")
    if df_genes_not_match.gene.tolist():
        lst_filter_gene_names += df_genes_not_match.gene.tolist()
        warnings_warn("No match between gene name from input file and retrieved gene name.\n{}".format(df_genes_not_match))
    return(lst_filter_gene_names)


def filter_genotypes(df_genotypes, df_phenotypes, lst_filter_gene_names=None, lst_filter_rs_id=None):
    lst_filter_ids_gt = (
        df_genotypes.loc[~df_genotypes.rs_id.str.startswith("rs", na=False)]
        .genotype_id
        .unique()
        .tolist()
    )  # no SVs
    if lst_filter_gene_names:
        lst_filter_ids_gt += df_genotypes.loc[df_genotypes.gene.isin(lst_filter_gene_names)].genotype_id.unique().tolist()
    if lst_filter_rs_id:
        lst_filter_ids_gt += df_genotypes.loc[df_genotypes.rs_id.isin(lst_filter_rs_id)].genotype_id.unique().tolist()
    lst_filter_ids_pt = df_phenotypes.loc[df_phenotypes.phenotype_name.isna()].genotype_id.unique().tolist()
    lst_filter_ids = lst_filter_ids_gt + lst_filter_ids_pt
    if len:
        print("Following genotype IDs (n={len}) are removed. {ids}".format(len=len(lst_filter_ids), ids=lst_filter_ids))
    df_filter_phenotypes = df_phenotypes[~df_phenotypes.genotype_id.isin(lst_filter_ids)]
    df_filter_genotypes = df_genotypes[~df_genotypes.genotype_id.isin(lst_filter_ids)]
    if df_filter_phenotypes.empty and df_filter_genotypes.empty:
        raise Exception("All genotypes are removed.")
    return(df_filter_phenotypes, df_filter_genotypes)


def write_bedfile(df_data, output_prefix, output_path):
    df_data_sub = df_data[["chrom", "start", "end", "name"]].sort_values(by="chrom", key=lambda x: np_argsort(
        index_natsorted(zip(df_data.chrom, df_data.start))))

    df_data_sub.to_csv("{path}{output_prefix}{ext}".format(
        path=str(output_path) + "/",
        output_prefix=output_prefix,
        ext=".bed",
        ), sep="\t", index=False, header=False)


def get_forward_orientation(variant_genotype):
    dict_complementary = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', '-': '.'}
    genotype = ""
    for base in variant_genotype:
        complementary = base.replace(base, dict_complementary.get(base, base))
        genotype = ''.join([genotype, complementary])
    return(genotype)


def get_indel_notation_with_flanking_base(record_ref, record_alt, variant_genotype):
    flanking_base_orientation = {
        "homozygote_ref": record_ref + "/" + record_ref,
        "heterozygote_ref_alt": record_ref + "/" + record_alt,
        "homozygote_alt": record_alt + "/" + record_alt,
    }
    ref = record_ref.ljust(len(record_alt), ".")
    alt = record_alt.ljust(len(record_ref), ".")
    for ref_base, alt_base in zip(ref, alt):
        if ref_base != alt_base:
            translate_notation = {
                ref_base + "/" + ref_base: "homozygote_ref",
                ref_base + "/" + alt_base: "heterozygote_ref_alt",
                alt_base + "/" + alt_base: "homozygote_alt",
            }
            if variant_genotype not in translate_notation:
                raise KeyError("Variant genotype type not expected.")
            genotype_simplified = translate_notation[variant_genotype]
            found_notation = flanking_base_orientation[genotype_simplified]
            break
        # print("From {} to {} for {} {} using ref {} and alt {}".format(
        #         variant_genotype,
        #         found_notation, record.rs_id, genotype_simplified, record_ref, record_alt
        #     )
        # )
    return(found_notation)


def generate_yaml_dict(df_phenotypes, df_metadata_variant_gt):
    yaml_dict = df_phenotypes.to_dict('index')
    for i, record in df_metadata_variant_gt.iterrows():
        if "snp" not in yaml_dict[record.genotype_id].keys():
            yaml_dict[record.genotype_id]["snp"] = list()
        # If genotype is based on a gene on the reverse strand, genotype is translated to the forward strand.
        if int(record.strand) == -1:
            record.variant_genotype = get_forward_orientation(variant_genotype=record.variant_genotype)
        if (
            (len(record.alt) > 1 and "/" not in record.alt)  # insertion
            or (len(record.ref) > 1 and "/" not in record.ref)  # deletion
        ):
            record.variant_genotype = get_indel_notation_with_flanking_base(
                record_ref=record.ref,
                record_alt=record.alt,
                variant_genotype=record.variant_genotype
            )
        yaml_dict[record.genotype_id]["snp"].append(
            record[["chrom", "start", "end", "rs_id", "name", "variant_genotype", "ref", "alt", "strand"]].to_dict()
        )
    return(yaml_dict)


def write_yaml(yaml_dict, output_prefix, output_path):
    with open("{path}/{prefix}.yaml".format(path=str(output_path) + "/", prefix=output_prefix), "w") as file:
        yaml.dump(yaml_dict, file, default_flow_style=False)


def compare_files(old, new):
    if isinstance(old, dict) and isinstance(new, dict):
        deepdiff_out = DeepDiff(old, new, ignore_order=True, verbose_level=2, report_repetition=True).pretty()
        if not deepdiff_out:
            print("Files are the same.")
        else:
            print("Differences between dictionaries:\n{}".format(deepdiff_out))
    elif isinstance(old, IOBase) and isinstance(new, IOBase):
        diff = unified_diff(old.readlines(), new.readlines(), n=0)
        delta = ''.join(x for x in diff)
        if not delta:
            print("Files are the same.")
        else:
            print(delta)
    else:
        warnings_warn("Comparison is not supported for this data type.")


def main(prev_bed_file, prev_yaml_file, output_path, output_prefix, config, translation_table):
    if not output_prefix:
        output_prefix = pathlib.Path(translation_table).stem.lower()
    df_phenotypes, df_genotypes = morph_translation_file(
        csv_file=translation_table,
        dict_rename_cols=config.getjsonloads("dict_rename_tf_cols")
    )
    df_ens_metadata = get_data_ensembl(
        df_translation=df_genotypes.loc[
            df_genotypes['rs_id'].str.startswith('rs', na=False),
            ['rs_id', 'gene_and_rs_id', 'gene']
        ].drop_duplicates(), 	# Select variant rows and relevant metadata columns.
        ensembl_url=config.get("ensembl_url"),
        species=config.get("species")
    )
    lst_filter_gene_names = get_invalid_genes_and_warn(
        df_ens_metadata=df_ens_metadata,
        genes_regex=config.getjsonloads("filter_gene_regex")
    )
    df_phenotypes, df_genotypes = filter_genotypes(
        df_genotypes=df_genotypes,
        df_phenotypes=df_phenotypes,
        lst_filter_gene_names=lst_filter_gene_names,
        lst_filter_rs_id=config.getjsonloads("lst_filter_rs_id", None)
    )
    write_bedfile(df_data=df_ens_metadata, output_path=str(output_path) + "/", output_prefix=output_prefix)
    if prev_bed_file:
        compare_files(old=prev_bed_file, new=open(output_path + output_prefix + ".bed"))

    # df_sv = (
    #     df_translation[
    #         ~df_translation['rs_id'].str.startswith('rs', na=False)
    #         & df_translation['genotype_id'].str.match("[0-9]")
    #     ]
    # )

    df_metadata_variant_gt = (
        pd
        .merge(
            df_ens_metadata,
            df_genotypes,
            how="left",
            left_on=["name", "rs_id", "gene"],
            right_on=["gene_and_rs_id", "rs_id", "gene"]
        )
        .dropna(axis=0, subset=["genotype_id"])
    )
    output_yaml = generate_yaml_dict(df_phenotypes=df_phenotypes, df_metadata_variant_gt=df_metadata_variant_gt)
    write_yaml(yaml_dict=output_yaml, output_prefix=output_prefix, output_path=output_path)
    if prev_yaml_file:
        prev_yaml = yaml.load(prev_yaml_file, Loader=yaml.FullLoader)
        compare_files(old=prev_yaml, new=output_yaml)


if __name__ == '__main__':
    args = parse_arguments_and_check(args_in=sys.argv[1:])
    config_section = read_config_section_and_check(section=args.config_section, config_file=args.config_file)
    main(
        prev_bed_file=args.bed,
        prev_yaml_file=args.yaml,
        output_path=args.output_path,
        output_prefix=args.output_prefix,
        config=config_section,
        translation_table=args.translation_table,
    )
