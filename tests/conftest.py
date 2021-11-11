import os
import pysam
import pytest
from pytest_reqs import check_requirements
import shutil
# @pytest.fixture(scope="session")
# def tmpdir(tmp_path_factory):
#     test_tmp_path = tmp_path_factory.mktemp("data")
#     return test_tmp_path
