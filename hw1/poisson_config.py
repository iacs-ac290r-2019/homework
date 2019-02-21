"""
Michael S. Emanuel
Wed Feb 20 21:35:04 2019
"""

import yaml

with open("example_1.yml", 'r') as cfg:
    try:
        print(yaml.load(cfg))
    except yaml.YAMLError as err:
        print('error in YAML.')
        print(err)
        
