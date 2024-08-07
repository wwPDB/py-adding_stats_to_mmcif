#!/usr/bin/env python
import argparse
import logging
import pprint

from .mmcif_dictionary_handling import SoftwareClassification
from .xml_parsing import parse_xml

logger = logging.getLogger()


# AL: not obvious decision, keeping as it was.
#all = 'Overall'
_all = ''

# AL:  corrections/additions in: pdbx_chi_squared, number_unique_obs,
# meanI_over_sigI_obs, pdbx_chi_squared, percent_possible_all
stats_remap = {
    'reflns': {'pos_list': ['Overall'],
               'cif_to_xml':
                   {'d_resolution_low': 'ResolutionLow',
                    'd_resolution_high': 'ResolutionHigh',
                    'pdbx_Rmerge_I_obs': 'Rmerge' + _all,
                    'pdbx_Rrim_I_all': 'Rmeas' + _all,
                    'pdbx_Rpim_I_all': 'Rpim' + _all,
                    'number_obs': 'NumberReflections',
                    'pdbx_netI_over_sigmaI': 'MeanIoverSD',
                    'pdbx_CC_half': 'CChalf',
                    'pdbx_chi_squared': 'MeanChiSq',
                    'percent_possible_obs': 'Completeness',
                    'pdbx_redundancy': 'Multiplicity'}},
    'reflns_shell': {'pos_list': ['Inner', 'Outer'],
                     'cif_to_xml':
                         {'d_res_low': 'ResolutionLow',
                          'd_res_high': 'ResolutionHigh',
                          'Rmerge_I_obs': 'Rmerge' + _all,
                          'pdbx_Rrim_I_all': 'Rmeas' + _all,
                          'pdbx_Rpim_I_all': 'Rpim' + _all,
                          'number_unique_obs': 'NumberReflections',
                          'meanI_over_sigI_obs': 'MeanIoverSD',
                          'pdbx_CC_half': 'CChalf',
                          'pdbx_chi_squared': 'MeanChiSq',
                          'percent_possible_all': 'Completeness',
                          'pdbx_redundancy': 'Multiplicity'}}
}

extra_cif_items = {'pdbx_ordinal': ''}

table_keys = {'deposition': 'reflns'}


def split_on_separator(separator, value_to_split):
    value_list = list()
    logging.debug('separator: "%s"' % separator)
    if separator:
        if separator == ' ':
            value_list = value_to_split.split()
        else:
            value_list = value_to_split.split(separator)
        logging.debug(value_list)
    return value_list


class aimlessReport:
    def __init__(self, xml_file):
        """
        :param xml_file: input aimless XML file
        """
        self.xml_file = xml_file
        self.tree = None
        self.root = None
        self.stats_dict = dict()
        self.cif_item = None
        self.cif_cat = None
        self.number_of_values = 0

    def parse_xml(self):
        """
            checks input file is XML file, parses it.
            prevents parsing twice by checking if self.tree already exists
        :return: True if a parsed aimless XML file, False if not
        """
        if not self.tree:
            self.root = parse_xml(xml_file=self.xml_file)

        if self.root is not None:
            if self.root.tag == 'AIMLESS_PIPE' or self.root.tag == 'AIMLESS':
                logging.debug('is an aimless xml file')
                return True
        return False

    def add_data_to_stats_dict(self, instance, value):
        self.stats_dict.setdefault(self.cif_cat, {}).setdefault(self.cif_item, [''] * self.number_of_values)[
            instance] = value

    def get_data(self):
        """
        :return: statistics dictionary from Dataset XML tag.
        """
        datasetresultnodes = self.root.findall(".//Result/Dataset")
        data_set_counter = 0
        for datasetresultnode in datasetresultnodes:
            data_set_counter += 1
            for self.cif_cat in stats_remap:

                location_list = stats_remap[self.cif_cat]['pos_list']
                self.number_of_values = len(location_list)

                for instance, location in enumerate(location_list):
                    for self.cif_item in stats_remap[self.cif_cat]['cif_to_xml']:
                        logging.debug(self.cif_item)
                        xml_item = stats_remap[self.cif_cat]['cif_to_xml'][self.cif_item]

                        xml_node = datasetresultnode.find(xml_item)
                        xml_item_for_location = xml_node.find(location)

                        logging.debug(xml_item_for_location)
                        xml_value = xml_item_for_location.text.strip()
                        logging.debug(xml_value)

                        self.add_data_to_stats_dict(instance=instance, value=xml_value)

                    for self.cif_item in extra_cif_items:
                        if extra_cif_items[self.cif_item]:
                            value = extra_cif_items[self.cif_item]
                        else:
                            value = instance + 1
                        self.add_data_to_stats_dict(instance=instance, value=str(value))

        return self.stats_dict

    def get_data_from_table(self):
        """
        :return: dictionary of statistics from CCP4 tables.
        """
        ccp4tables = self.root.findall(".//CCP4Table")
        logging.debug(ccp4tables)
        for table in ccp4tables:
            logging.debug(table.attrib)
            if 'id' in table.attrib:
                if table.attrib['id'] in table_keys:
                    # need to set cif category based on the table name.
                    self.cif_cat = table_keys[table.attrib['id']]
                    headers = table.find('headers')
                    header_list = split_on_separator(separator=headers.attrib['separator'].text,
                                                     value_to_split=headers.text)
                    data = table.find('data').text
                    logging.debug(data)
                    data_lines = data.strip().split('\n')
                    logging.debug(data_lines)
                    number_of_data_values = len(data_lines)
                    logging.debug('number of data items: %s' % number_of_data_values)
                    for instance, d in enumerate(data_lines):
                        d = split_on_separator(separator=headers.attrib['separator'].text, value_to_split=d)
                        for header_pos, self.cif_item in enumerate(header_list):
                            logging.debug('%s - position %s' % (self.cif_item, header_pos))
                            value = d[header_pos]
                            self.add_data_to_stats_dict(instance=instance, value=value)

        return self.stats_dict

    def return_data(self):
        """
        master function which returns data from aimless XML file.
        :return:
        """
        is_aimless_file = self.parse_xml()
        if is_aimless_file:
            self.get_data_from_table()
            if not self.stats_dict:
                self.get_data()

        return self.stats_dict

    def get_aimlesss_version(self):
        """
        :return: aimless version from the XML file
        """
        version = None
        is_aimless_file = self.parse_xml()
        if is_aimless_file:
            header = self.root.findall(".//AIMLESS")
            for head in header:
                if 'version' in head.attrib:
                    version = head.attrib['version']
        return version

    def get_aimless_version_dict(self):
        version = self.get_aimlesss_version()
        software_row = SoftwareClassification().get_software_row(software_name='aimless', version=version)

        return software_row


if __name__ == '__main__':  # pragma: no cover
    parser = argparse.ArgumentParser()
    parser.add_argument('--xml_file', help='input xml file', type=str, required=True)
    parser.add_argument('-d', '--debug', help='debugging', action='store_const', dest='loglevel', const=logging.DEBUG,
                        default=logging.INFO)

    args = parser.parse_args()
    logger.setLevel(args.loglevel)

    xml_file = args.xml_file
    ar = aimlessReport(xml_file=xml_file)
    xml_data = ar.return_data()
    pprint.pprint(xml_data)
