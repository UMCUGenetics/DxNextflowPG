#! venv/bin/python
# standard libraries alphabetic order of main package.
import configparser
import csv
import shutil
import pandas as pd

# third party libraries alphabetic order of main package.
import requests
import pytest

# local libraries alphabetic order of main package.
import assets.create_workflow_reference_files as create_ref


# TODO implement these tests.
"""
rsID and/or location linked to multiple genes -> not unique for single genotype gene.
variant with multiple alternative alleles ensembl
variant with multiple alternative alleles translation table
"""


@pytest.fixture(scope="module")
def create_test_files(tmp_path_factory):
    test_path = str(tmp_path_factory.mktemp("tests")) + "/"
    shutil.rmtree(test_path)
    shutil.copytree("tests/test_data/", test_path)
    return test_path

class TestCreateRefsParser():
    def test_parser_optional_config_section(self):
        parser = create_ref.parse_arguments_and_check(
            args_in=["--config_section", "fake"]
        )
        assert parser
    
    def test_parser_optional_bed(self):
        parser = create_ref.parse_arguments_and_check(
            args_in=["--bed", "references/sites_of_interest.bed"]
        )
        assert parser

    def test_parser_optional_output_path(self):
        parser = create_ref.parse_arguments_and_check(
            args_in=["--output_path", "references/"]
        )
        assert parser

    def test_parser_optional_output_prefix(self):
        parser = create_ref.parse_arguments_and_check(
            args_in=["--output_prefix", "fake_prefix"]
        )
        assert parser

    def test_parser_optional_yaml(self):
        parser = create_ref.parse_arguments_and_check(
            args_in=["--yaml", "references/sites_of_interest.yaml"]
        )
        assert parser

    def test_parser_optional_bed_not_exist(self):
        with pytest.raises(FileNotFoundError):
            parser = create_ref.parse_arguments_and_check(
                args_in=["--bed", "non_existing_file.bed"]
            )

    def test_parser_optional_output_path_not_exist(self):
        with pytest.raises(FileNotFoundError):
            parser = create_ref.parse_arguments_and_check(
                args_in=["--output_path", "non_existing_path/"]
            )

    def test_parser_optional_yaml_not_exist(self):
        with pytest.raises(FileNotFoundError):
            parser = create_ref.parse_arguments_and_check(
                args_in=["--yaml", "non_existing_file.yaml"]
            )


class TestCreateRefsConfig():
    def test_config_file_not_exists(self):
        with pytest.raises(FileNotFoundError):
            config_section = create_ref.read_config_section_and_check(
                section="DEFAULT", config_file="non_existing_file.ini"
            )

    def test_config_section_not_exists(self, create_test_files):
        with pytest.raises(KeyError):
            config_section = create_ref.read_config_section_and_check(
                section="non_existing_section", config_file=create_test_files + "/ini_config/correct_config.ini"
            )

    def test_config_key_not_exists(self, create_test_files):
        with pytest.raises(KeyError):
            config_section = create_ref.read_config_section_and_check(
                section="DEFAULT", config_file=create_test_files + "/ini_config/missing_key.ini"
            )

    def test_config_file_empty(self, create_test_files):
        with pytest.raises(OSError) as empty_error:
            config_section = create_ref.read_config_section_and_check(
                section="DEFAULT", config_file=create_test_files + "/ini_config/empty.ini"
            )
        assert "empty" in str(empty_error.value)
    
    def test_config_translation_empty(self, create_test_files):
        with pytest.raises(OSError) as empty_error:
            config_section = create_ref.read_config_section_and_check(
                section="DEFAULT", config_file=create_test_files + "/ini_config/empty_translation_file.ini"
            )
        assert "empty" in str(empty_error.value)

    def test_config_translation_not_exists(self, create_test_files):
        with pytest.raises(FileNotFoundError):
            config_section = create_ref.read_config_section_and_check(
                section="DEFAULT", config_file=create_test_files + "/ini_config/non_existing_translation_file.ini"
            )


class TestCreateRefsTranslationFile():
    def test_tt_correct(self, create_test_files):
        tt_pheno, tt_geno = create_ref.morph_translation_file(csv_file=create_test_files + "/translation_input_files_extern/tt_correct.csv")
        assert True

    def test_tt_rename_columns_not_exists(self, create_test_files):
        tt_pheno, tt_geno = create_ref.morph_translation_file(
            csv_file=create_test_files + "/translation_input_files_extern/tt_correct.csv",
            dict_rename_cols={"non_existing_column": "fake_rename"}
        )
        assert True
    
    def test_tt_rename_columns(self, create_test_files):
        tt_pheno, tt_geno = create_ref.morph_translation_file(
            csv_file=create_test_files + "/translation_input_files_extern/tt_rename_columns.csv",
            dict_rename_cols={"variant_id": "gene_and_rs_id"}
        )
        assert True

    def test_tt_missing_gt_id(self, create_test_files):
        tt_pheno, tt_geno = create_ref.morph_translation_file(csv_file=create_test_files + "/translation_input_files_extern/tt_missing_genotype_id.csv")
        assert True

    def test_tt_missing_phenotypes(self, create_test_files):
        tt_pheno, tt_geno = create_ref.morph_translation_file(csv_file=create_test_files + "/translation_input_files_extern/tt_missing_phenotypes.csv")
        assert True

    def test_tt_missing_variant_gts(self, create_test_files):
        with pytest.raises(ValueError) as err_missing_variant_gt:
            tt_pheno, tt_geno = create_ref.morph_translation_file(
                csv_file=create_test_files + "/translation_input_files_extern/tt_missing_variant_genotypes.csv")
        assert "Variant_genotype value is required" in str(err_missing_variant_gt.value)

    def test_tt_variant_gt_format(self, create_test_files):
        tt_pheno, tt_geno = create_ref.morph_translation_file(csv_file=create_test_files + "/translation_input_files_extern/tt_variant_genotype_format.csv")
        assert True

    def test_tt_variant_gt_unexpected_sep(self, create_test_files):
        with pytest.raises(ValueError) as err_unexpected_sep:
            tt_pheno, tt_geno = create_ref.morph_translation_file(
                csv_file=create_test_files + "/translation_input_files_extern/tt_variant_genotype_unexpected_sep.csv"
            )
        assert "Expected single separator ':' or '/'" in str(err_unexpected_sep.value)
        assert "variant_genotype" in str(err_unexpected_sep.value)
    
    def test_tt_variant_gt_invalid_char(self, create_test_files):
        with pytest.raises(ValueError) as err_invalid_char:
            tt_pheno, tt_geno = create_ref.morph_translation_file(
                csv_file=create_test_files + "/translation_input_files_extern/tt_variant_genotype_invalid_char.csv"
            )
        assert "Invalid character in variant genotype. Supported:" in str(err_invalid_char.value)

    def test_tt_variant_id_unexpected_sep(self, create_test_files):
        with pytest.raises(ValueError) as err_unexpected_sep:
            tt_pheno, tt_geno = create_ref.morph_translation_file(
                csv_file=create_test_files + "/translation_input_files_extern/tt_variant_id_unexpected_sep.csv"
            )
        assert "Expected separator _" in str(err_unexpected_sep.value)
        assert "gene_and_rs_id" in str(err_unexpected_sep.value)

    def test_tt_variant_id_multi_sep(self, create_test_files):
        with pytest.raises(ValueError) as err_unexpected_sep:
            tt_pheno, tt_geno = create_ref.morph_translation_file(
                csv_file=create_test_files + "/translation_input_files_extern/tt_variant_id_multi_sep.csv")
        assert "Expected separator _" in str(err_unexpected_sep.value)
        assert "gene_and_rs_id" in str(err_unexpected_sep.value)

    def test_tt_missing_columns(self, create_test_files):
        with pytest.raises(ValueError) as err_missing_cols:
            tt_pheno, tt_geno = create_ref.morph_translation_file(csv_file=create_test_files + "/translation_input_files_extern/tt_columns_incomplete.csv")
        assert "Required columns are missing" in str(err_missing_cols.value)
    

class TestCreateRefsEnsembl():
    def test_ensembl(self):
        query_result = create_ref.get_ensembl_request_response(
            server="http://grch37.rest.ensembl.org", ext="/variation/homo_sapiens/", json={"ids": ["rs1799853"]}, method="post"
        )
        assert query_result
    
    def test_ensembl_typo(self):
        with pytest.raises(requests.exceptions.ConnectionError):
            query_result = create_ref.get_ensembl_request_response(
                server="http://grch37.rest.ensembl_typo.org", ext="/variation/homo_sapiens/", json={"ids": ["rs1799853"]}, method="post"
            )

    def test_rs_id_not_linked_to_gene(self):
        with pytest.warns(UserWarning, match="Gene not found"):
            create_ref.get_gene_metadata_ensembl(server="http://grch37.rest.ensembl.org",
                                                 ens_rs_id="rs1633021", location="6:29746868-29746869", species="homo_sapiens")

    def test_rs_id_not_linked_to_gene(self):
        with pytest.warns(UserWarning, match="Variation identifiers not found in Ensembl"):
            create_ref.get_variant_metadata_ensembl(server="http://grch37.rest.ensembl.org",
                                                    ids=["fakers"], species="homo_sapiens")


class TestCreateRefsInvalidGenes():
    def test_gene_with_variants_of_diff_chroms(self):
        df_metadata = pd.DataFrame.from_dict(
            {'row_1': ["fakegene", "fakegene", 1], 'row_2': ["fakegene", "fakegene", 2], },
            orient='index',
            columns=['gene', 'retrieved_gene', 'chrom']
        )
        with pytest.warns(UserWarning, match="At least one gene is linked to variants from different chromosomes."):
            create_ref.get_invalid_genes_and_warn(df_ens_metadata=df_metadata, genes_regex=None)

    def test_gene_no_match_retrieved(self):
        df_metadata = pd.DataFrame.from_dict(
            {'row_1': ["fakegene", "fakegene2", 1]},
            orient='index',
            columns=['gene', 'retrieved_gene', 'chrom']
        )
        with pytest.warns(UserWarning, match="No match between gene name from input file and retrieved gene name"):
            create_ref.get_invalid_genes_and_warn(df_ens_metadata=df_metadata, genes_regex=None)


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
            df_filter_pt, df_filter_gt = create_ref.filter_genotypes(df_genotypes=df_genotypes, df_phenotypes=df_phenotypes)
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


# TODO implement these tests in TestCreateRefsVariantGenotypes
class TestCreateRefsVariantGenotypeNotation():
    def test_forward_oriantation(self):
        dict_forward_reverse_genotypes = {
            "A/A": "T/T", # hom ref
            "A/T": "T/A", # ref alt
            "T/A": "A/T", # alt ref
            "C/G": "G/C", # remaining bases
            "./A": "./T", # insertion
            "A/.": "T/.", # deletion
        }
        for gt_forward, gt_reverse in dict_forward_reverse_genotypes.items():
            retrieved_gt = create_ref.get_forward_orientation(variant_genotype=gt_reverse)
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
            notation = create_ref.get_indel_notation_with_flanking_base(
                record_ref="AA", record_alt="A", variant_genotype="./A"
            )
        assert "Variant genotype type not expected." in str(wrong_notation.value)


# TODO implement these tests in TestCreateRefsCompareFiles
"""
compare files, not supported data type.
compare files, dict yaml, only order different dict
compare files, dict yaml, different variant 
compare files, dict yaml, removed genotype
compare files, dict yaml, removed random field.
compare files, bed, removed rs site.
compare files, bed, synonym rs id and thus name.
"""
class TestCreateRefsCompareFiles():
    pass


# TODO implement these tests in TestCreateRefsYaml
"""
yaml write, file exists already. - -> ERROR, param 'force'?
yaml write, file exists already, use force - -> succeed.
variant genotypes 1 bp on reverse strand are changed to forward notation.
- REF/REF
- REF/ALT
- ALT/REF
- ALT/ALT
variant genotypes > 1bp on reverse strand are changed to forward notation.
- REF/REF
- REF/ALT
- ALT/REF
- ALT/ALT
variants genotypes indels notation standardized.
- REF/REF
- REF/ALT
- ALT/REF
- ALT/ALT
"""
class TestCreateRefsYaml():
    pass


# TODO implement these tests in TestCreateRefsBed
"""
write bedfile to non existing file - -> succeed
write bedfile, file exists already - -> ERROR, param 'force'?
"""
class TestCreateRefsBed():
    pass


# TODO implement these tests in TestCreateRefsOutput
"""
output file starts with output_prefix.
output file saved in output path.
"""
class TestCreateRefsOutput():
    pass

