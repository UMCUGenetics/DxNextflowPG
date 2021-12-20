#! venv/bin/python
# standard libraries alphabetic order of main package.
import pathlib
import pandas as pd
from pathlib import Path
import shutil

# third party libraries alphabetic order of main package.
import pytest
import requests
import yaml

# local libraries alphabetic order of main package.
import assets.create_workflow_reference_files as create_ref


@pytest.fixture(scope="module")
def create_test_files(tmp_path_factory):
    test_path = str(tmp_path_factory.mktemp("tests")) + "/"
    shutil.rmtree(test_path)
    shutil.copytree("tests/test_data/", test_path)
    return test_path


@pytest.fixture(scope="module")
def get_test_yaml(create_test_files):
    with open(create_test_files + "/translation_yaml/sites_of_interest.yaml") as yaml_file:
        prev_yaml = yaml.load(yaml_file, Loader=yaml.FullLoader)
    return prev_yaml


@pytest.fixture(scope="module")
def get_tmp_output_path(tmp_path_factory):
    output_path = str(tmp_path_factory.mktemp("output")) + "/"
    open(str(output_path) + "existing_output.bed", "a").close()
    open(str(output_path) + "existing_output.yaml", "a").close()
    return output_path


class TestCreateRefsParser():
    def test_parser_translation_table(self, create_test_files):
        parser = create_ref.parse_arguments_and_check(
            args_in=[create_test_files + "/translation_input_files_extern/tt_correct.csv"]
        )
        assert parser.translation_table

    def test_parser_translation_table_not_exist(self):
        with pytest.raises(FileNotFoundError) as not_exist_error:
            create_ref.parse_arguments_and_check(
                args_in=["non_existing_tt.csv"]
            )
        assert "No such file or directory" in str(not_exist_error.value)

    def test_parser_translation_table_empty(self, create_test_files):
        with pytest.raises(OSError) as empty_error:
            create_ref.parse_arguments_and_check(
                args_in=[create_test_files + "/translation_input_files_extern/tt_empty.csv"]
            )
        assert "empty" in str(empty_error.value)

    def test_parser_optional_config_file(self, create_test_files):
        parser = create_ref.parse_arguments_and_check(
            args_in=[
                "--config_file", "./assets/create_workflow_reference_files.ini",
                create_test_files + "/translation_input_files_extern/tt_correct.csv"
            ]
        )
        assert parser

    def test_parser_optional_config_file_not_exists(self, create_test_files):
        with pytest.raises(FileNotFoundError) as not_exist_error:
            create_ref.parse_arguments_and_check(
                args_in=[
                    "--config_file", "non_existing_file.ini",
                    create_test_files + "/translation_input_files_extern/tt_correct.csv"
                ]
            )
        assert "No such file or directory" in str(not_exist_error.value)

    def test_parser_optional_config_file_empty(self, create_test_files):
        with pytest.raises(OSError) as empty_error:
            create_ref.parse_arguments_and_check(
                args_in=[
                    "--config_file", create_test_files + "/ini_config/empty.ini",
                    create_test_files + "/translation_input_files_extern/tt_correct.csv"
                ]
            )
        assert "empty" in str(empty_error.value)

    def test_parser_optional_config_section(self, create_test_files):
        parser = create_ref.parse_arguments_and_check(
            args_in=[
                "--config_section", "fake",
                create_test_files + "/translation_input_files_extern/tt_correct.csv"
            ]
        )
        assert parser

    def test_parser_optional_bed(self, create_test_files):
        parser = create_ref.parse_arguments_and_check(
            args_in=[
                "--bed", "references/sites_of_interest_GRCh38.bed",
                create_test_files + "/translation_input_files_extern/tt_correct.csv"
            ]
        )
        assert parser

    def test_parser_optional_output_path_absolute(self, create_test_files):
        parser = create_ref.parse_arguments_and_check(
            args_in=[
                "--output_path", str(pathlib.Path.cwd()),
                create_test_files + "/translation_input_files_extern/tt_correct.csv"
            ]
        )
        assert parser.output_path

    def test_parser_optional_output_path_relative(self, create_test_files):
        with pytest.raises(OSError) as absolute_path_error:
            create_ref.parse_arguments_and_check(
                args_in=[
                    "--output_path", "references/",
                    create_test_files + "/translation_input_files_extern/tt_correct.csv"
                ]
            )
        assert "absolute" in str(absolute_path_error.value)

    def test_parser_optional_output_prefix(self, create_test_files):
        parser = create_ref.parse_arguments_and_check(
            args_in=[
                "--output_prefix", "fake_prefix",
                create_test_files + "/translation_input_files_extern/tt_correct.csv"
            ]
        )
        assert parser

    def test_parser_optional_yaml(self, create_test_files):
        parser = create_ref.parse_arguments_and_check(
            args_in=[
                "--yaml", "references/sites_of_interest_GRCh38.yaml",
                create_test_files + "/translation_input_files_extern/tt_correct.csv"
            ]
        )
        assert parser

    def test_parser_optional_bed_not_exist(self, create_test_files):
        with pytest.raises(FileNotFoundError) as not_exist_error:
            create_ref.parse_arguments_and_check(
                args_in=[
                    "--bed", "non_existing_file.bed",
                    create_test_files + "/translation_input_files_extern/tt_correct.csv"
                ]
            )
        assert "No such file or directory" in str(not_exist_error.value)

    def test_parser_optional_output_path_not_exist(self, create_test_files):
        with pytest.raises(OSError) as path_not_exist_error:
            create_ref.parse_arguments_and_check(
                args_in=[
                    "--output_path", "non_existing_path/",
                    create_test_files + "/translation_input_files_extern/tt_correct.csv"
                ]
            )
        print(str(path_not_exist_error.value))
        assert "No such file or directory" in str(path_not_exist_error.value)
        assert "non_existing_path" in str(path_not_exist_error.value)

    def test_parser_optional_yaml_not_exist(self, create_test_files):
        with pytest.raises(FileNotFoundError) as not_exist_error:
            create_ref.parse_arguments_and_check(
                args_in=[
                    "--yaml", "non_existing_file.yaml",
                    create_test_files + "/translation_input_files_extern/tt_correct.csv"
                ]
            )
        assert "No such file or directory" in str(not_exist_error.value)


class TestCreateRefsConfig():
    def test_config_file(self):
        config_section = create_ref.read_config_section_and_check(
            section="DEFAULT", config_file=open("./assets/create_workflow_reference_files.ini")
        )
        assert config_section

    def test_config_section_not_exists(self, create_test_files):
        with pytest.raises(KeyError) as keyerror:
            create_ref.read_config_section_and_check(
                section="non_existing_section", config_file=open(create_test_files + "/ini_config/correct_config.ini")
            )
        assert "non_existing_section" in str(keyerror.value)

    def test_config_key_not_exists(self, create_test_files):
        with pytest.raises(KeyError) as keyerror:
            create_ref.read_config_section_and_check(
                section="DEFAULT", config_file=open(create_test_files + "/ini_config/missing_key.ini")
            )
        print(str(keyerror.value))
        assert "Required key ensembl_url not in config file" in str(keyerror.value)


class TestCreateRefsTranslationFile():
    def test_tt_correct(self, create_test_files):
        tt_pheno, tt_geno = create_ref.morph_translation_file(
            csv_file=create_test_files + "/translation_input_files_extern/tt_correct.csv"
        )
        assert not tt_pheno.empty
        assert not tt_geno.empty

    def test_tt_correct_with_sections(self, create_test_files):
        tt_pheno, tt_geno = create_ref.morph_translation_file(
            csv_file=create_test_files + "/translation_input_files_extern/tt_correct_with_sections.csv"
        )
        assert not tt_pheno.empty
        assert not tt_geno.empty

    def test_tt_genotype_with_sv(self, create_test_files):
        tt_pheno, tt_geno = create_ref.morph_translation_file(
            csv_file=create_test_files + "/translation_input_files_extern/tt_genotype_with_sv.csv"
        )
        assert not tt_pheno.empty
        assert not tt_geno.empty

    def test_tt_rename_columns_not_exists(self, create_test_files):
        tt_pheno, tt_geno = create_ref.morph_translation_file(
            csv_file=create_test_files + "/translation_input_files_extern/tt_correct.csv",
            dict_rename_cols={"non_existing_column": "fake_rename"}
        )
        assert not tt_pheno.empty
        assert not tt_geno.empty

    def test_tt_rename_columns(self, create_test_files):
        tt_pheno, tt_geno = create_ref.morph_translation_file(
            csv_file=create_test_files + "/translation_input_files_extern/tt_rename_columns.csv",
            dict_rename_cols={"variant_id": "gene_and_rs_id"}
        )
        assert not tt_pheno.empty
        assert not tt_geno.empty

    def test_tt_missing_gt_id(self, create_test_files):
        tt_pheno, tt_geno = create_ref.morph_translation_file(
            csv_file=create_test_files + "/translation_input_files_extern/tt_missing_genotype_id.csv"
        )
        assert not tt_pheno.empty
        assert not tt_geno.empty

    def test_tt_missing_phenotypes(self, create_test_files):
        tt_pheno, tt_geno = create_ref.morph_translation_file(
            csv_file=create_test_files + "/translation_input_files_extern/tt_missing_phenotypes.csv"
        )
        assert not tt_pheno.empty
        assert not tt_geno.empty

    def test_tt_missing_variant_gts(self, create_test_files):
        with pytest.raises(ValueError) as err_missing_variant_gt:
            tt_pheno, tt_geno = create_ref.morph_translation_file(
                csv_file=create_test_files + "/translation_input_files_extern/tt_missing_variant_genotypes.csv"
            )
        assert "Variant_genotype value is required" in str(err_missing_variant_gt.value)

    def test_tt_variant_gt_format(self, create_test_files):
        tt_pheno, tt_geno = create_ref.morph_translation_file(
            csv_file=create_test_files + "/translation_input_files_extern/tt_variant_genotype_format.csv"
        )
        assert not tt_pheno.empty
        assert not tt_geno.empty

    def test_tt_variant_gt_unexpected_sep(self, create_test_files):
        with pytest.raises(ValueError) as err_unexpected_sep:
            create_ref.morph_translation_file(
                csv_file=create_test_files + "/translation_input_files_extern/tt_variant_genotype_unexpected_sep.csv"
            )
        assert "Expected single separator ':' or '/'" in str(err_unexpected_sep.value)
        assert "variant_genotype" in str(err_unexpected_sep.value)

    def test_tt_variant_gt_invalid_char(self, create_test_files):
        with pytest.raises(ValueError) as err_invalid_char:
            create_ref.morph_translation_file(
                csv_file=create_test_files + "/translation_input_files_extern/tt_variant_genotype_invalid_char.csv"
            )
        assert "Invalid character in variant genotype. Supported:" in str(err_invalid_char.value)

    def test_tt_variant_id_unexpected_sep(self, create_test_files):
        with pytest.raises(ValueError) as err_unexpected_sep:
            create_ref.morph_translation_file(
                csv_file=create_test_files + "/translation_input_files_extern/tt_variant_id_unexpected_sep.csv"
            )
        assert "Expected separator _" in str(err_unexpected_sep.value)
        assert "gene_and_rs_id" in str(err_unexpected_sep.value)

    def test_tt_variant_id_multi_sep(self, create_test_files):
        with pytest.raises(ValueError) as err_unexpected_sep:
            create_ref.morph_translation_file(
                csv_file=create_test_files + "/translation_input_files_extern/tt_variant_id_multi_sep.csv"
            )
        assert "Expected separator _" in str(err_unexpected_sep.value)
        assert "gene_and_rs_id" in str(err_unexpected_sep.value)

    def test_tt_missing_columns(self, create_test_files):
        with pytest.raises(ValueError) as err_missing_cols:
            create_ref.morph_translation_file(
                csv_file=create_test_files + "/translation_input_files_extern/tt_columns_incomplete.csv"
            )
        assert "Required columns are missing" in str(err_missing_cols.value)


class TestCreateRefsEnsembl():
    def test_ensembl(self):
        query_result = create_ref.get_ensembl_request_response(
            server="http://rest.ensembl.org", ext="/variation/homo_sapiens/", json={"ids": ["rs1799853"]}, method="post"
        )
        assert query_result

    def test_ensembl_typo(self):
        with pytest.raises(requests.exceptions.ConnectionError):
            create_ref.get_ensembl_request_response(
                server="http://rest.ensembl_typo.org",
                ext="/variation/homo_sapiens/",
                json={"ids": ["rs1799853"]},
                method="post"
            )

    def test_rs_id_not_linked_to_gene(self):
        with pytest.warns(UserWarning, match="Gene not found"):
            create_ref.get_gene_metadata_ensembl(
                server="http://rest.ensembl.org",
                ens_rs_id="rs1633021",
                location="chr6:29746868-29746869",
                species="homo_sapiens"
            )

    def test_rs_id_not_exists(self):
        with pytest.warns(UserWarning, match="Variation identifiers not found in Ensembl"):
            create_ref.get_variant_metadata_ensembl(
                server="http://rest.ensembl.org", ids=["fakers"], species="homo_sapiens"
            )


class TestCreateRefsInvalidGenes():
    def test_gene_with_variants_of_diff_chroms(self):
        df_metadata = pd.DataFrame.from_dict(
            {'row_1': ["fakegene", "fakegene", 'chr1'], 'row_2': ["fakegene", "fakegene", 'chr2'], },
            orient='index',
            columns=['gene', 'retrieved_gene', 'chrom']
        )
        with pytest.warns(UserWarning, match="At least one gene is linked to variants from different chromosomes."):
            create_ref.get_invalid_genes_and_warn(df_ens_metadata=df_metadata, genes_regex=None)

    def test_gene_no_match_retrieved(self):
        df_metadata = pd.DataFrame.from_dict(
            {'row_1': ["fakegene", "fakegene2", "chr1"]},
            orient='index',
            columns=['gene', 'retrieved_gene', 'chrom']
        )
        with pytest.warns(UserWarning, match="No match between gene name from input file and retrieved gene name"):
            create_ref.get_invalid_genes_and_warn(df_ens_metadata=df_metadata, genes_regex=None)


# TODO: remove sites with multiple alternatives (based on Ensembl) aka not duploid
class TestCreateRefsFilterGenotypes():
    def test_filter_genotypes_succes(self, create_test_files):
        df_phenotypes, df_genotypes = create_ref.morph_translation_file(
            csv_file=create_test_files + "/translation_input_files_extern/tt_correct.csv"
        )
        df_filter_pt, df_filter_gt = create_ref.filter_genotypes(df_genotypes=df_genotypes, df_phenotypes=df_phenotypes)
        assert not df_filter_pt.empty
        assert not df_filter_gt.empty

    def test_filter_genotypes_sv(self, create_test_files):
        # I wonder if this dependency should be removed?
        df_phenotypes, df_genotypes = create_ref.morph_translation_file(
            csv_file=create_test_files + "/translation_input_files_extern/tt_genotype_with_sv.csv"
        )
        with pytest.raises(Exception) as genotypes_removed:
            df_filter_pt, df_filter_gt = create_ref.filter_genotypes(
                df_genotypes=df_genotypes, df_phenotypes=df_phenotypes
            )
        assert "All genotypes are removed." in str(genotypes_removed.value)

    def test_filter_genotypes_remove_gene(self, create_test_files):
        # I wonder if this dependency should be removed?
        df_phenotypes, df_genotypes = create_ref.morph_translation_file(
            csv_file=create_test_files + "/translation_input_files_extern/tt_genotype_rs_multi.csv"
        )
        with pytest.raises(Exception) as genotypes_removed:
            df_filter_pt, df_filter_gt = create_ref.filter_genotypes(
                df_genotypes=df_genotypes, df_phenotypes=df_phenotypes, lst_filter_gene_names=["fakegene"]
            )
        assert "All genotypes are removed." in str(genotypes_removed.value)

    def test_filter_genotypes_remove_rs_id(self, create_test_files):
        # I wonder if this dependency should be removed?
        df_phenotypes, df_genotypes = create_ref.morph_translation_file(
            csv_file=create_test_files + "/translation_input_files_extern/tt_genotype_rs_multi.csv"
        )
        with pytest.raises(Exception) as genotypes_removed:
            df_filter_pt, df_filter_gt = create_ref.filter_genotypes(
                df_genotypes=df_genotypes, df_phenotypes=df_phenotypes, lst_filter_rs_id=["rs0000001"]
            )
        assert "All genotypes are removed." in str(genotypes_removed.value)


class TestCreateRefsVariantGenotypeNotation():
    def test_complementary_sequence(self):
        dict_forward_reverse_genotypes = {
            "A/A": "T/T",  # hom ref
            "A/T": "T/A",  # ref alt
            "T/A": "A/T",  # alt ref
            "C/G": "G/C",  # remaining bases
            "./A": "./T",  # insertion
            "A/.": "T/.",  # deletion
        }
        for gt_forward, gt_reverse in dict_forward_reverse_genotypes.items():
            retrieved_gt = create_ref.get_complementary_sequence(variant_genotype=gt_reverse)
            assert retrieved_gt == gt_forward

    def test_indel_notation_ins_ref(self):
        notation = create_ref.get_indel_notation_with_flanking_base(
            record_ref="A", record_alt="AA", variant_genotype="./."
        )
        assert notation == "A/A"

    def test_indel_notation_ins_het(self):
        notation = create_ref.get_indel_notation_with_flanking_base(
            record_ref="A", record_alt="AA", variant_genotype="./A"
        )
        assert notation == "A/AA"

    def test_indel_notation_ins_hom_alt(self):
        notation = create_ref.get_indel_notation_with_flanking_base(
            record_ref="A", record_alt="AA", variant_genotype="A/A"
        )
        assert notation == "AA/AA"

    def test_indel_notation_del_ref(self):
        notation = create_ref.get_indel_notation_with_flanking_base(
            record_ref="AA", record_alt="A", variant_genotype="A/A"
        )
        assert notation == "AA/AA"

    def test_indel_notation_del_het(self):
        notation = create_ref.get_indel_notation_with_flanking_base(
            record_ref="AA", record_alt="A", variant_genotype="A/."
        )
        assert notation == "AA/A"

    def test_indel_notation_del_alt(self):
        notation = create_ref.get_indel_notation_with_flanking_base(
            record_ref="AA", record_alt="A", variant_genotype="./."
        )
        assert notation == "A/A"

    def test_indel_notation_wrong_genotype(self):
        with pytest.raises(KeyError) as wrong_notation:
            create_ref.get_indel_notation_with_flanking_base(
                record_ref="AA", record_alt="A", variant_genotype="./A"
            )
        assert "Variant genotype type not expected." in str(wrong_notation.value)


class TestCreateRefsCompareFiles():
    def test_cf_yaml_same_files(self, get_test_yaml, capsys):
        create_ref.compare_dict(old=get_test_yaml, new=get_test_yaml)
        captured = capsys.readouterr()
        assert captured.out == "Files are the same.\n"

    def test_cf_yaml_removed_snp_site(self, create_test_files, get_test_yaml, capsys):
        with open(create_test_files + "/translation_yaml/removed_snp_site.yaml") as yaml_file:
            dict_new = yaml.load(yaml_file, Loader=yaml.FullLoader)
        create_ref.compare_dict(old=get_test_yaml, new=dict_new)
        captured = capsys.readouterr()
        assert "Differences between dictionaries:" in captured.out
        assert "root['576']['snp'][1] removed" in captured.out
        assert "root['577']['snp'][1]" in captured.out
        assert "root['578']['snp'][1]" in captured.out

    def test_cf_yaml_removed_snp_field_name(self, create_test_files, get_test_yaml, capsys):
        with open(create_test_files + "/translation_yaml/removed_snp_field_name.yaml") as yaml_file:
            dict_new = yaml.load(yaml_file, Loader=yaml.FullLoader)
        create_ref.compare_dict(old=get_test_yaml, new=dict_new)
        captured = capsys.readouterr()
        assert "Differences between dictionaries:" in captured.out
        for gt in ['576', '577', '578']:
            assert "root['{gt}']['snp'][0] removed".format(gt=gt) in captured.out
            assert "root['{gt}']['snp'][1] removed".format(gt=gt) in captured.out

    def test_cf_yaml_removed_genotype(self, create_test_files, get_test_yaml, capsys):
        with open(create_test_files + "/translation_yaml/removed_genotype.yaml") as yaml_file:
            dict_new = yaml.load(yaml_file, Loader=yaml.FullLoader)
        create_ref.compare_dict(old=get_test_yaml, new=dict_new)
        captured = capsys.readouterr()
        assert "Differences between dictionaries:" in captured.out
        assert "root['576'] removed" in captured.out

    def test_cf_yaml_diff_phenotype_metadata(self, create_test_files, get_test_yaml, capsys):
        with open(create_test_files + "/translation_yaml/diff_phenotype_metadata.yaml") as yaml_file:
            dict_new = yaml.load(yaml_file, Loader=yaml.FullLoader)
        create_ref.compare_dict(old=get_test_yaml, new=dict_new)
        captured = capsys.readouterr()
        assert "Differences between dictionaries:" in captured.out
        assert "Value of root['576']['phenotype_id'] changed from 534.0 to 1534.0" in captured.out

    def test_cf_yaml_diff_order_snp(self, create_test_files, get_test_yaml, capsys):
        with open(create_test_files + "/translation_yaml/diff_order_snp.yaml") as yaml_file:
            dict_new = yaml.load(yaml_file, Loader=yaml.FullLoader)
        create_ref.compare_dict(old=get_test_yaml, new=dict_new)
        captured = capsys.readouterr()
        assert captured.out == "Files are the same.\n"

    def test_cf_yaml_diff_order_genotypes(self, create_test_files, get_test_yaml, capsys):
        with open(create_test_files + "/translation_yaml/diff_order_genotypes.yaml") as yaml_file:
            dict_new = yaml.load(yaml_file, Loader=yaml.FullLoader)
        create_ref.compare_dict(old=get_test_yaml, new=dict_new)
        captured = capsys.readouterr()
        assert captured.out == "Files are the same.\n"

    def test_cf_yaml_diff_genotype_metadata(self, create_test_files, get_test_yaml, capsys):
        with open(create_test_files + "/translation_yaml/diff_genotype_metadata.yaml") as yaml_file:
            dict_new = yaml.load(yaml_file, Loader=yaml.FullLoader)
        create_ref.compare_dict(old=get_test_yaml, new=dict_new)
        captured = capsys.readouterr()
        assert "Differences between dictionaries:" in captured.out
        for gt in ['576', '577', '578']:
            for idx in [0, 1]:
                assert "root['{gt}']['snp'][{i}] added to iterable.".format(gt=gt, i=idx) in captured.out
                assert "root['{gt}']['snp'][{i}] removed from iterable.".format(gt=gt, i=idx) in captured.out
            assert "Value of root['{gt}']['gene'] changed from".format(gt=gt) in captured.out
            assert "Value of root['{gt}']['gene_genotype'] changed from".format(gt=gt), captured.out
            assert "TPMT2" in captured.out

    def test_cf_yaml_added_snp_site(self, create_test_files, get_test_yaml, capsys):
        with open(create_test_files + "/translation_yaml/added_snp_site.yaml") as yaml_file:
            dict_new = yaml.load(yaml_file, Loader=yaml.FullLoader)
        create_ref.compare_dict(old=get_test_yaml, new=dict_new)
        captured = capsys.readouterr()
        assert "Differences between dictionaries:" in captured.out
        for gt in ['576', '577', '578']:
            assert "root['{gt}']['snp'][2] added to iterable.".format(gt=gt) in captured.out

    def test_cf_yaml_added_genotype(self, create_test_files, get_test_yaml, capsys):
        with open(create_test_files + "/translation_yaml/added_genotype.yaml") as yaml_file:
            dict_new = yaml.load(yaml_file, Loader=yaml.FullLoader)
        create_ref.compare_dict(old=get_test_yaml, new=dict_new)
        captured = capsys.readouterr()
        assert "Differences between dictionaries:" in captured.out
        assert "root['579'] added" in captured.out

    def test_cf_bed_same_file(self, create_test_files, capsys):
        create_ref.compare_tab_files(
            old=open(create_test_files + "/bed_files/sites_of_interest.bed"),
            new=open(create_test_files + "/bed_files/sites_of_interest.bed")
        )
        captured = capsys.readouterr()
        assert captured.out == "Files are the same.\n"

    def test_cf_bed_removed_site(self, create_test_files, capsys):
        create_ref.compare_tab_files(
            old=open(create_test_files + "/bed_files/sites_of_interest.bed"),
            new=open(create_test_files + "/bed_files/removed_site.bed")
        )
        captured = capsys.readouterr()
        assert "-chr6\t18130917\t18130918\tTPMT_rs1142345" in captured.out

    def test_cf_bed_diff_pos(self, create_test_files, capsys):
        create_ref.compare_tab_files(
            old=open(create_test_files + "/bed_files/sites_of_interest.bed"),
            new=open(create_test_files + "/bed_files/diff_pos.bed")
        )
        captured = capsys.readouterr()
        assert "-chr6\t18139227\t18139228\tTPMT_rs1800460" in captured.out
        assert "+chr6\t18139200\t18139201\tTPMT_rs1800460" in captured.out

    def test_cf_bed_diff_name(self, create_test_files, capsys):
        create_ref.compare_tab_files(
            old=open(create_test_files + "/bed_files/sites_of_interest.bed"),
            new=open(create_test_files + "/bed_files/diff_name.bed")
        )
        captured = capsys.readouterr()
        assert "-chr6\t18139227\t18139228\tTPMT_rs1800460" in captured.out
        assert "+chr6\t18139227\t18139228\tTPMT_rs1800000" in captured.out

    def test_cf_bed_added_site(self, create_test_files, capsys):
        create_ref.compare_tab_files(
            old=open(create_test_files + "/bed_files/sites_of_interest.bed"),
            new=open(create_test_files + "/bed_files/added_site.bed")
        )
        captured = capsys.readouterr()
        assert "+chr6\t18143954\t18143955\tTPMT_rs1800462" in captured.out


class TestCreateRefsYaml():
    def test_generate_yaml_rev_strand(self):
        df_data_rev = pd.DataFrame.from_dict(
            {
                'row_1': ['1', "chr6", 1000, 1001, "fakeRS1", "name", "T/T", "A", "G", -1],
                'row_2': ['1', "chr6", 1001, 1002, "fakeRS2", "name", "T/C", "A", "G", -1],
                'row_3': ['1', "chr6", 1002, 1003, "fakeRS3", "name", "C/T", "A", "G", -1],
                'row_4': ['1', "chr6", 1003, 1004, "fakeRS4", "name", "C/C", "A", "G", -1],
            },
            orient='index',
            columns=["genotype_id", "chrom", "start", "end", "rs_id", "name", "variant_genotype", "ref", "alt", "strand"]
        )
        df_pheno = pd.DataFrame.from_dict(
            {'1': ['1', "fakegene", "wildtype/wildtype", "fakegene:wildtype/wildtype", 1, "POOR METABOLIZER"]},
            orient='index',
            columns=["genotype_id", "gene", "genotype_realname", "gene_genotype", "phenotype_id", "phenotype_name"]
        )
        yaml_dict = create_ref.generate_translation_dict(df_phenotypes=df_pheno, df_metadata_variant_gt=df_data_rev)
        assert '1' in yaml_dict
        assert yaml_dict['1']['snp'][0]['variant_genotype'] == 'A/A'
        assert yaml_dict['1']['snp'][1]['variant_genotype'] == 'A/G'
        assert yaml_dict['1']['snp'][2]['variant_genotype'] == 'G/A'
        assert yaml_dict['1']['snp'][3]['variant_genotype'] == 'G/G'

    def test_generate_yaml_rev_strand_indel(self):
        df_data_rev = pd.DataFrame.from_dict(
            {
                'row_1': ['1', "chr6", 1000, 1001, "fakeRS1", "name", "./.", "A", "AA", -1],
                'row_2': ['1', "chr6", 1001, 1002, "fakeRS2", "name", "./T", "A", "AA", -1],
                'row_3': ['1', "chr6", 1002, 1003, "fakeRS3", "name", "T/T", "A", "AA", -1],
                'row_4': ['2', "chr6", 1000, 1001, "fakeRS1", "name", "./.", "AA", "A", -1],
                'row_5': ['2', "chr6", 1002, 1003, "fakeRS3", "name", "T/.", "AA", "A", -1],
                'row_6': ['2', "chr6", 1003, 1004, "fakeRS4", "name", "T/T", "AA", "A", -1],
            },
            orient='index',
            columns=["genotype_id", "chrom", "start", "end", "rs_id", "name", "variant_genotype", "ref", "alt", "strand"]
        )
        df_pheno = pd.DataFrame.from_dict(
            {'1': ['1', "fake_gene_ins", "wildtype/wildtype", "fake_gene_ins:wildtype/wildtype", 1, "POOR METABOLIZER"],
             '2': ['2', "fake_gene_del", "wildtype/wildtype", "fake_gene_del:wildtype/wildtype", 1, "POOR METABOLIZER"]},
            orient='index',
            columns=["genotype_id", "gene", "genotype_realname", "gene_genotype", "phenotype_id", "phenotype_name"]
        )
        yaml_dict = create_ref.generate_translation_dict(df_phenotypes=df_pheno, df_metadata_variant_gt=df_data_rev)
        assert '1' in yaml_dict
        assert '2' in yaml_dict
        assert yaml_dict['1']['snp'][0]['variant_genotype'] == 'A/A'
        assert yaml_dict['1']['snp'][1]['variant_genotype'] == 'A/AA'
        assert yaml_dict['1']['snp'][2]['variant_genotype'] == 'AA/AA'
        assert yaml_dict['2']['snp'][0]['variant_genotype'] == 'A/A'
        assert yaml_dict['2']['snp'][1]['variant_genotype'] == 'AA/A'
        assert yaml_dict['2']['snp'][2]['variant_genotype'] == 'AA/AA'


class TestCreateRefsWriteBed():
    def test_write_bed(self, get_tmp_output_path):
        df_data = pd.DataFrame.from_dict(
            {'row_1': ["chr1", 1000, 1001, "test"]},
            orient='index',
            columns=['chrom', 'start', 'end', 'name']
        )
        create_ref.write_bedfile(df_data=df_data, output_prefix="test", output_path=get_tmp_output_path)
        assert Path(get_tmp_output_path + "test.bed").is_file()

    def test_write_bed_existing_file(self, get_tmp_output_path):
        df_data = pd.DataFrame.from_dict(
            {'row_1': ["chr1", 1000, 1001, "test"]},
            orient='index',
            columns=['chrom', 'start', 'end', 'name']
        )
        create_ref.write_bedfile(df_data=df_data, output_prefix="existing_output", output_path=get_tmp_output_path)
        assert Path(get_tmp_output_path + "existing_output.bed").is_file()


class TestCreateRefsWriteYaml():
    def test_write_yaml(self, get_tmp_output_path):
        create_ref.write_yaml(yaml_dict={"key1": "value1"}, output_prefix="test", output_path=get_tmp_output_path)
        assert Path(get_tmp_output_path + "test.yaml").is_file()

    def test_write_yaml_existing_file(self, get_tmp_output_path):
        create_ref.write_yaml(yaml_dict={"key1": "value1"}, output_prefix="existing_output", output_path=get_tmp_output_path)
        assert Path(get_tmp_output_path + "existing_output.yaml").is_file()
