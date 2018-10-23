#!/usr/bin/env python
import argparse
import drsip

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='DRSIP: Distance Restraints- and Cyclic Symmetry-Imposed Packing.\n\nPredicts the quaternary structure of Homo-oligomeric Transmembrane Proteins (HoTPs) assuming the complex is Cn-symmetric.')
    parser.add_argument(
        'static_pdb', help='Static monomer PDB file')
    parser.add_argument(
        'mobile_pdb', help='Mobile monomer PDB file')
    parser.add_argument('zdock_output_file',
                        help='ZDOCK 3.0.2\'s output file')
    parser.add_argument('-d', dest='dist_rest_file',
                        help='Distance restraints file', required=False, default='')
    parser.add_argument('-o', dest='output_file',
                        help='Output the top N results file. Default: \'DRSIP_results.xlsx\'. See -N', required=False, default='DRSIP_results.xlsx')
    parser.add_argument('-p', dest='poses_folder',
                        help='Folder to write the top N poses. Default: Current folder. See -N', required=False, default='')
    parser.add_argument('-N', dest='num_poses',
                        help='Set "N", the number of top poses to write out. Default: 10', required=False, type=int, default=10)
    parser.add_argument('--soluble', dest='soluble', help='Use DRSIP\'s soluble protein docking protocol. When calling DRSIP with --soluble, --dssp and --trans-helix are not required and have no effect. Optional.',
                        required=False, action='store_true', default=False)
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('--dssp', dest='dssp_path',
                       help='Path to the DSSP program. Required if --trans-helix is not provided. Cannot set --dssp and --trans-helix at the same time. Set one or the other', default='')
    group.add_argument('--trans-helix', dest='trans_helix', help='Comma separated resids of transmembrane helices.\
    Ex. Helix 1: resid 5-20 , Helix 2: resid 23-40, Helix 3: resid 45-60 then TRANS-HELIX:"5-20, 23-40, 45-60". These can be obtained from the Orientations of Proteins in Membranes (OPM) database or PPM server. Cannot set --dssp and --trans-helix at the same time. Set one or the other', default='')
    parser.add_argument('--save-all', dest='save_file', help='For advanced users only. Saves all the input files and filtering data to a compressed msgpack file. Optional.',
                        required=False, default='')
    args = parser.parse_args()

    if args.soluble:

        if args.dist_rest_file == '':
            parser.error(
                "--soluble requires a valid distance restraints file set with -d")

    else:

        if len(args.dssp_path) == 0 and len(args.trans_helix) == 0:
            parser.error(
                "The transmembrane protocol requires the path of the DSSP program (--dssp) or the list of transmembrane helices (--trans-helix)")

    drsip_obj = drsip.DR_SIP(args.static_pdb, args.mobile_pdb,
                             args.zdock_output_file, args.dist_rest_file, args.soluble)

    if not args.soluble:

        if len(args.trans_helix) > 0:
            helical_elements_sel_str = ['resid %s' %
                                        x.strip() for x in args.trans_helix.split()]

            drsip_obj.run(helical_elements_sel_str=helical_elements_sel_str)

        else:
            drsip_obj.run(dssp_path=args.dssp_path)

    else:
        drsip_obj.run()

    if len(args.save_file) != 0:
        drsip_obj.save_data(args.save_file)

    drsip_obj.final_results.to_excel(args.output_file)
    drsip_obj.write_topN_poses(args.poses_folder, num_poses=args.num_poses)