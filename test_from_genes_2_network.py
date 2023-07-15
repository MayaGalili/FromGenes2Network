import unittest
from from_genes_2_network import FromGenes2Networks


class TestFromGenesToNetwork(unittest.TestCase):

    # Each method in this class that starts with 'test_' will be run as a test case
    def test_whole_run(self):
        net_generator = FromGenes2Networks("sequence.gb")
        net_generator.run()