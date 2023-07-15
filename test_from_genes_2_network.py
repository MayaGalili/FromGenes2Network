import os.path
import unittest
from from_genes_2_network import FromGenes2Networks


class TestFromGenesToNetwork(unittest.TestCase):

    def test_whole_run(self):
        script_dir = os.path.dirname(os.path.realpath(__file__))

        net_generator = FromGenes2Networks(os.path.join(script_dir, "sequence.gb"), script_dir)
        net_generator.run()
        self.assertTrue(os.path.isfile(os.path.join(script_dir, 'GenesAlignmentNetwork.png')))