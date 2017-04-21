import unittest
from reaction import react_str

class TestMols(unittest.TestCase):
    def test_basic(self):
        self.assertEqual(react_str("CCO"), ["CC=O","CC(=O)O"])

    def test_two(self):
        self.assertEqual(react_str("CCC"), [])

    def test_three(self):
        self.assertEqual(react_str("CCCS"), [])


if __name__ == '__main__':
    unittest.main()