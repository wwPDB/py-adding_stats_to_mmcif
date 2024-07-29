#!/usr/bin/env python
import argparse
import logging
import os

from .aimless_xml_parser import aimlessReport
from .cif_handling import mmcifHandling

logger = logging.getLogger()
FORMAT = "%(filename)s - %(funcName)s - %(message)s"
logging.basicConfig(format=FORMAT)


def get_xml_data(xml_file):
    # get data from aimless XML file
    ar = aimlessReport(xml_file=xml_file)
    xml_data = ar.return_data()
    software_row = ar.get_aimless_version_dict()

    return xml_data, software_row

def fix_resolution_limits(pc):
    # fix discrepancy between low resolution limits
    b = pc.cif_handling.cifObj[0]
    r1 = b.find_values('_refine.ls_d_res_low')
    r2 = b.find_values('_refine_ls_shell.d_res_low')
    ai = b.find_values('_reflns.d_resolution_low')[0]
    b.find_values('_reflns_shell.d_res_low')[0] = ai
    if float(r1[0]) > float(ai):
      r1[0] = ai
    if float(r2[0]) > float(ai):
      r2[0] = ai
    r1 = b.find_values('_refine.ls_d_res_high')
    r2 = b.find_values('_refine_ls_shell.d_res_high')
    ai = b.find_values('_reflns.d_resolution_high')[-1]
    b.find_values('_reflns_shell.d_res_high')[-1] = ai
    if float(r1[-1]) < float(ai):
      r1[-1] = ai
    if float(r2[-1]) < float(ai):
      r2[-1] = ai

def fix_resolution_cross_val(pc):
    b = pc.cif_handling.cifObj[0]
    c = b.get_mmcif_category('_refine')
    for k in c:
      c[k] = c[k][0]
    c['pdbx_ls_cross_valid_method'] = 'FREE R-VALUE'
    b.set_pairs('_refine.', c, False)

def run_process(xml_file, input_cif, output_cif):
    xml_data, software_row = get_xml_data(xml_file=xml_file)
    if xml_data:
        # if there is data from the XML file then add this to the mmCIF file
        if os.path.exists(input_cif):
            pc = mmcifHandling()
            pc.parse_mmcif(fileName=input_cif)
            # add aimless data to the mmCIF file
            ok = pc.addToCif(data_dictionary=xml_data)
            fix_resolution_limits(pc)
            fix_resolution_cross_val(pc)
            # update the software list in the mmCIF file to add aimless
            software_cat = pc.addValuesToCategory(category='software', item_value_dictionary=software_row,
                                                  ordinal_item='pdbx_ordinal')
            ok = pc.addToCif(data_dictionary=software_cat)
            # add exptl data
            pc.addExptlToCif()
            # write out the resulting mmCIF file.
            pc.writeCif(fileName=output_cif)
            if os.path.exists(output_cif):
                return True
    return False


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_mmcif', help='output mmcif file', type=str, required=True)
    parser.add_argument('--input_mmcif', help='input mmcif file', type=str, required=True)
    parser.add_argument('--xml_file', help='input xml file', type=str, required=True)
    parser.add_argument('-d', '--debug', help='debugging', action='store_const', dest='loglevel', const=logging.DEBUG,
                        default=logging.INFO)

    args = parser.parse_args()

    logger.setLevel(args.loglevel)

    run_process(xml_file=args.xml_file, input_cif=args.input_mmcif, output_cif=args.output_mmcif)
