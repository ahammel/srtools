from srtools import stats
from srtools.test import TEST_FOLDER
import os

def test_summary_statistics():
    stats.print_summary_statistics(TEST_FOLDER + "/test_data/speed_test.sam", output_file="test_summary.txt")
    test_file = open("test_summary.txt")
    known_file = open(TEST_FOLDER + "/test_data/david_summary.txt")
    test_lines = test_file.readlines()
    known_lines = known_file.readlines()
    for line in test_lines:
        assert line in known_lines
    for line in known_lines:
        assert line in test_lines
    test_file.close()
    known_file.close()
    os.remove("test_summary.txt")


