import os
import pysam
import pytest
from pytest_reqs import check_requirements
import shutil
# @pytest.fixture(scope="session")
# def tmpdir(tmp_path_factory):
#     test_tmp_path = tmp_path_factory.mktemp("data")
#     return test_tmp_path


@pytest.fixture(scope="session", autouse=True)
def setup_and_get_test_path(tmp_path_factory):
    test_tmp_path = tmp_path_factory.mktemp("data")

    # create empty files
    open(str(test_tmp_path) + "/empty.yaml", "a").close()
    open(str(test_tmp_path) + "/empty.json", "a").close()
    open(str(test_tmp_path) + "/empty.vcf.gz", "a").close()
    open(str(test_tmp_path) + "/empty.vcf.gz.tbi", "a").close()
    for vcf_file in ["./tests/test_data/no_records.vcf", "./tests/test_data/fake_000000000000_R00C00.vcf"]:
        shutil.copy(vcf_file, test_tmp_path)
        basename = os.path.basename(vcf_file)
        tmp_vcf_file = str(test_tmp_path) + "/" + basename
        pysam.tabix_compress(tmp_vcf_file, tmp_vcf_file + ".gz")
        pysam.tabix_index(tmp_vcf_file + ".gz", preset="vcf")
    shutil.copy(str(test_tmp_path) + "/" + "fake_000000000000_R00C00.vcf",
                str(test_tmp_path) + "/" + "fake_000000000001_R01C01.vcf")
    print(os.listdir(test_tmp_path))
    return str(test_tmp_path) + "/"

