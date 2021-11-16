#! venv/bin/python
# import statements, alphabetic order of main package.
import argparse
from configparser import ConfigParser
from difflib import unified_diff
from errno import ENOENT as errno_ENOENT
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

def parse_arguments_and_check(args_in):
	parser = argparse.ArgumentParser(
		description="Migrate translation table to a bed file with sites of interest " +
		"and a yaml file with genotype - phenotype translation.")
	parser.add_argument("-s", "--config_section", type=str, default="DEFAULT", help="The INI-config section used.")
	parser.add_argument("-b", "--bed", type=str, required=False,
	                    help="Previous bed file. If provided, differences between previous and generated bed file are shown.")
	parser.add_argument("-o", "--output_path", type=str, required=False, default=pathlib.Path().resolve(),
	                    help="Filepath where output is placed. (default: %(default)s)")
	parser.add_argument("-p", "--output_prefix", type=str, required=False,
						help="Output prefix to use as filename in output. Default is filename prefix of translation table.")
	parser.add_argument("-y", "--yaml", type=str, required=False,
						help="Previous yaml file. If provided, differences between previous and generated yaml file are shown.")
	args = parser.parse_args(args_in)

	if args.bed:
		check_file(file=args.bed)
	if args.yaml:
		check_file(file=args.yaml)
	if args.output_path:
		check_file(file=args.output_path)
	return(args)


def check_file(file):
	if not pathlib.Path(file).is_file() and not pathlib.Path(file).is_dir():
		raise FileNotFoundError(errno_ENOENT, os_strerror(errno_ENOENT), file)
	elif not pathlib.Path(file).stat().st_size:
		raise OSError("File is empty.")


def read_config_section_and_check(section, config_file="./assets/create_workflow_reference_files.ini"):
	config_parser = ConfigParser(converters={"jsonloads": json.loads})
	check_file(file=config_file)
	config_parser.read(config_file)
	config_section = config_parser[section]
	required_keys = ["ensembl_url", "species", "translation_table"]
	for req_key in required_keys:
		if req_key not in config_section:
			raise KeyError("Required key {} not in config file.".format(req_key))
	file_keys = ["translation_table"]
	for file_key in file_keys:
		check_file(config_section.get(file_key))
	return(config_section)


def morph_input_file(csv_file, dict_rename_cols):
	df_translation = pd.read_csv(csv_file).rename(columns=dict_rename_cols)
	df_translation[['gene', 'rs_id']] = df_translation.gene_and_rs_id.str.split("_", expand=True,)
	df_translation['variant_genotype'] = df_translation['variant_genotype'].str.replace(':', "/")
	# split translation table into phenotype and genotype metadata.
	df_phenotypes = (
		df_translation[["genotype_id", "gene_genotype", "gene", "genotype_realname", "phenotype_id", "phenotype_name"]]
		.dropna(axis=0, subset=["gene", "phenotype_name"])
		.set_index("genotype_id", drop=False)
	)
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


def get_gene_metadata_ensembl(server, ens_rs_id, location, species):
	client_out = get_ensembl_request_response(
		server=server,
		ext="/overlap/region/{species}/{loc}?feature=gene;logic_name=ensembl_havana_gene_homo_sapiens_37".format(
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
	get_variant_metadata_ensembl(server=ensembl_url, ids=df_translation.rs_id.tolist(), species=species)
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
				'chrom': map_info["seq_region_name"],
				'start': map_info["start"]-1, # transform to zero based.
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
	lst_filter_ids_gt = df_genotypes.loc[~df_genotypes.rs_id.str.startswith("rs", na=False)].genotype_id.unique().tolist()  # no SVs
	if lst_filter_gene_names:
		lst_filter_ids_gt += df_genotypes.loc[df_genotypes.gene.isin(lst_filter_gene_names)].genotype_id.unique().tolist()
	if lst_filter_rs_id:
		lst_filter_ids_gt += df_genotypes.loc[df_genotypes.rs_id.isin(lst_filter_rs_id)].genotype_id.unique().tolist()
	lst_filter_ids_pt = df_phenotypes.loc[df_phenotypes.phenotype_name.isna()].genotype_id.unique().tolist()
	lst_filter_ids = lst_filter_ids_gt + lst_filter_ids_pt
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
		path=output_path,
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
		"heterozygote": record_ref + "/" + record_alt,
		"homozygote_alt": record_alt + "/" + record_alt,
	}
	ref = record_ref.ljust(len(record_alt), ".")
	alt = record_alt.ljust(len(record_ref), ".")
	for ref_base, alt_base in zip(ref, alt):
		if ref_base != alt_base:
			translate_notation = {
				ref_base + "/" + ref_base: "homozygote_ref",
				ref_base + "/" + alt_base: "heterozygote",
				alt_base + "/" + alt_base: "homozygote_alt",
			}
			if variant_genotype not in translate_notation:
				raise KeyError("Variant genotype type not expected.")
			genotype_simplified = translate_notation[variant_genotype]
			found_notation = flanking_base_orientation[genotype_simplified]
		# print("From {} to {} for {} {} using ref {} and alt {}".format(variant_genotype,
        #                                                          found_notation, record.rs_id, genotype_simplified, record_ref, record_alt))
	return(found_notation)


def generate_yaml_dict(df_phenotypes, df_metadata_variant_gt, output_path, output_prefix):
	yaml_dict = df_phenotypes.to_dict('index')
	for i, record in df_metadata_variant_gt.iterrows():
		if "snp" not in yaml_dict[record.genotype_id].keys():
			yaml_dict[record.genotype_id]["snp"] = list()
		# If genotype is based on a gene on the reverse strand, genotype is translated to the forward strand.
		if int(record.strand) == -1:
			record.variant_genotype = get_forward_orientation(variant_genotype=record.variant_genotype)
		# if indel, standardize notation.
		if len(record.alt) > 1 and "/" not in record.alt and "." in record.variant_genotype:
			record.variant_genotype = get_indel_notation_with_flanking_base(
				record_ref=record.ref, 
				record_alt=record.alt, 
				variant_genotype=record.variant_genotype
			)
		yaml_dict[record.genotype_id]["snp"].append(
			record[["chrom", "start", "end", "rs_id", "name", "variant_genotype", "ref", "alt", "strand"]].to_dict()
			)
	with open(output_path + output_prefix + ".yaml", "w") as file:
		documents = yaml.dump(yaml_dict, file, default_flow_style=False)
	return(yaml_dict)


def compare_files(old, new):
	if isinstance(old, dict) and isinstance(new, dict):
		print("Differences between dictionaries:\n{}".format(DeepDiff(old, new, ignore_order=True).pretty()))
	elif (
		isinstance(old, str) 
		and pathlib.PurePath(old).suffix == ".bed"
        and isinstance(new, str) 
		and pathlib.PurePath(new).suffix == ".bed"
		):
		with open(old, 'r') as prev:
			with open(new, 'r') as current:
				diff = unified_diff(prev.readlines(), current.readlines(), fromfile='previous bed', tofile='current bed', n=0)
				print("Differences between bed files:")
				for line in diff:
					print(line)
	else:
		print("Comparison is not supported for this data type.")


def main(prev_bed_file, prev_yaml_file, output_path, output_prefix, config):
	check_file(file=config.get("translation_table"))
	if not args.output_prefix:
		args.output_prefix = pathlib.Path(config.get("translation_table")).stem.lower()
	df_phenotypes, df_genotypes = morph_input_file(
		csv_file=config.get("translation_table"),
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
		lst_filter_rs_id=config.getjsonloads("lst_filter_rs_id")
	)
	write_bedfile(df_data=df_ens_metadata, output_path=output_path, output_prefix=output_prefix)
	compare_files(old=prev_bed_file, new=output_path + output_prefix + ".bed")

	# df_sv = df_translation[~df_translation['rs_id'].str.startswith('rs', na=False) & df_translation['genotype_id'].str.match("[0-9]")]

	df_metadata_variant_gt = (
            pd
            .merge(df_ens_metadata, df_genotypes, how="left", left_on="name", right_on="gene_and_rs_id")
            .dropna(axis=0, subset=["genotype_id"])
        )
	output_yaml = generate_yaml_dict(
        df_phenotypes=df_phenotypes, 
		df_genotypes=df_genotypes,
		df_metadata_variant_gt=df_metadata_variant_gt,
		output_path=output_path, 
		output_prefix=output_prefix
	)
	with open(prev_yaml_file) as yaml_file:
		prev_yaml = yaml.load(yaml_file, Loader=yaml.FullLoader)
	compare_files(old=prev_yaml, new=output_yaml)
	


if __name__ == '__main__':
	args = parse_arguments_and_check(args_in=sys.argv[1:])
	config_section = read_config_section_and_check(section=args.config_section, config_file="./assets/create_workflow_reference_files.ini")
	main(
		prev_bed_file=args.bed, 
		prev_yaml_file=args.yaml,
		output_path=args.output_path, 
		output_prefix=args.output_prefix,
		config=config_section,
	)

