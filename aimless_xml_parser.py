#!/usr/bin/env python
import xml.etree.ElementTree as ET
import logging
import pprint

logger = logging.getLogger()

stats_remap = {
    'reflns': {'pos_list': ['Overall'],
               'cif_to_xml':
                   {'d_resolution_low': 'ResolutionLow',
                    'd_resolution_high': 'ResolutionHigh',
                    'pdbx_Rmerge_I_obs': 'Rmerge',
                    'pdbx_Rrim_I_all': 'Rmeas',
                    'pdbx_Rpim_I_all': 'Rpim',
                    'number_measured_obs': 'NumberReflections',
                    'pdbx_netI_over_sigmaI': 'MeanIoverSD',
                    'pdbx_CC_half': 'CChalf',
                    'percent_possible_obs': 'Completeness',
                    'pdbx_redundancy': 'Multiplicity'}},
    'reflns_shell': {'pos_list': ['Inner', 'Outer'],
                     'cif_to_xml':
                         {'d_res_low': 'ResolutionLow',
                          'd_res_high': 'ResolutionHigh',
                          'Rmerge_I_obs': 'Rmerge',
                          'pdbx_Rrim_I_all': 'Rmeas',
                          'pdbx_Rpim_I_all': 'Rpim',
                          'number_measured_obs': 'NumberReflections',
                          'pdbx_netI_over_sigmaI': 'MeanIoverSD',
                          'pdbx_CC_half': 'CChalf',
                          'percent_possible_obs': 'Completeness',
                          'pdbx_redundancy': 'Multiplicity'}}
}

extra_cif_items = {'pdbx_diffrn_id': '1',
                   'pdbx_ordinal': ''}

class aimlessReport:
    def __init__(self, xml_file):
        self.xml_file = xml_file
        self.tree = None
        self.root = None
        self.stats_dict = dict()

    def parse_xml(self):
        self.tree = ET.parse(self.xml_file)
        self.root = self.tree.getroot()

    def get_data(self):
        datasetresultnodes = self.root.findall(".//Result/Dataset")
        for datasetresultnode in datasetresultnodes:
            for cif_cat in stats_remap:

                location_list = stats_remap[cif_cat]['pos_list']
                number_of_values = len(location_list)

                for instance, location in enumerate(location_list):
                    for cif_item in stats_remap[cif_cat]['cif_to_xml']:
                        logging.debug(cif_item)
                        xml_item = stats_remap[cif_cat]['cif_to_xml'][cif_item]

                        xml_node = datasetresultnode.find(xml_item)
                        xml_item_for_location = xml_node.find(location)

                        logging.debug(xml_item_for_location)
                        xml_value = xml_item_for_location.text.strip()
                        logging.debug(xml_value)

                        self.stats_dict.setdefault(cif_cat, {}).setdefault(cif_item, ['']*number_of_values)[instance] = xml_value
                    for cif_item in extra_cif_items:
                        if extra_cif_items[cif_item]:
                            value = extra_cif_items[cif_item]
                        else:
                            value = instance + 1
                        self.stats_dict.setdefault(cif_cat, {}).setdefault(cif_item, [''] * number_of_values)[instance] = str(value)

        return self.stats_dict

if __name__ == '__main__':
    logger.setLevel(logging.DEBUG)
    xml_file = 'test_data/gam-pipe.xml'
    ar = aimlessReport(xml_file=xml_file)
    ar.parse_xml()
    ar.get_data()
    pprint.pprint(ar.stats_dict)
