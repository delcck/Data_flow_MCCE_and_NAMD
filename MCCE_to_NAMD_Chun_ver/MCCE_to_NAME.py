from MCCE_to_NAMD_Chun_ver.base.data import Info
from pathlib import Path
import pandas as pd
from collections import defaultdict
import os

class MSAnaToNAMD(object):
    '''
    Task: Convert a microstat analysis output from MCCE to
    structural files ready for NAMD runs.
    '''

    def __init__(
            self,
            in_ms_csv: Path=None,
            in_pdb: Path=None,
            in_topologies: list=None):
        '''
        Task: Initialize the object

        Args:
            in_ms_csv:
                - the input path for the CSV that holds MCCE microstate analyses results for coming NAMD run.

            in_pdb:
                - the input PDB file for MCCE calculation.

            in_topologies:
                - a list of paths pointing to topology files to b sourced.
        '''
        self.in_ms_csv = in_ms_csv
        self.in_ms_df = pd.read_csv(in_ms_csv)
        self.Info = Info()
        self.in_pdb = in_pdb
        self.in_topologies = in_topologies


    def get_example(self):
        '''
        Task: Provide an example on the format for the input CSV.
        '''
        example_dict = {
            'GLUA0007_': [0.0, 0.0, 0.0, 0.0, 0.0],
            'HISA0015_': [1.0, 1.0, 1.0, 1.0, 1.0],
            'ASPA0018_': [0.0, 0.0, 0.0, -1.0, 0.0],
            'TYRA0020_': [0.0, 0.0, 0.0, 0.0, 0.0],
            'CTRA0129_': [0.0, 0.0, 0.0, 0.0, 0.0],
            'Count': [1129301.0, 190004.0, 53465.0, 36175.0, 36096.0],
            'Occupancy': [0.753, 0.127, 0.036, 0.024, 0.024],
            'Sum_crg_protein': [19.0, 18.0, 18.0, 18.0, 18.0],
        }
        self.example = pd.DataFrame(example_dict)
        print(self.example.head())


    def parse_info_from_ms_ana(self, number_of_chosen_cases=None):
        '''
        Task: Parse information from the input MMCE calculations and convert them
        into actionable lines for users and psfgen in VMD.

        Args:
            number_of_chosen_cases
                - The number of cases to be chosen from MCCE microstate analyses results.
        '''
        #--1st find residues involved, residues types involved, and their charge states present
        self.__get_chain_and_residue_types__()
        self.__get_per_type_charge_states__()
        #--2nd select a subset of data as being instructed
        if not number_of_chosen_cases or number_of_chosen_cases >= len(self.in_ms_df):
            self.sel_ms_df = self.in_ms_df.copy(deep=True)
            number_of_chosen_cases = len(self.sel_ms_df)
        elif number_of_chosen_cases > 0:
            self.sel_ms_df = self.in_ms_df.head(number_of_chosen_cases).copy(deep=True)
        else:
            print('The "number_of_chosen_cases" must be set to be larger than 0.')
        pass
        #--3rd collect messages for each case present in self.sel_ms_df
        actions_dict = defaultdict(lambda: [])
        for ind in range(0, number_of_chosen_cases, 1):
            actions = self.__get_message_and_psfgen_lines_for_a_row__(
                self.sel_ms_df.loc[[ind]]
            )
            actions_dict[f'MCCE_case_{str(ind)}'] = actions
        #--4th group actions by chain id for psfgen processing
        self.actions_by_chains = self.__sort_actions_by_chains__(actions_dict)
        print('Information from the input MCCE calculations has been parsed and can be accessed by the variable "self.actions_by_chains".')


    def print_action_messages(self):
        for case_key in self.actions_by_chains.keys():
            print(f'Actions to be taken for {case_key} are .......')
            for chain_key in self.actions_by_chains[case_key].keys():
                print(f'Actions for chain {chain_key} of {case_key}:')
                for ind, line in enumerate(self.actions_by_chains[case_key][chain_key]):
                    print(f'Action {ind}: {line}')
            print('''
Next case ....
''')


    def write_psfgen_file(self, out_dir: str):
        '''
        Task: generate a compact TCL script for psfgen as to build the NAMD system from
        MCCE calculations.

        Args:
            out_dir:
                - the output directory to deposit psfgen files generated.
        '''
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        #--generate that for 1 case
        for case_key in self.actions_by_chains.keys():
            out_path = os.path.join(out_dir, f'{case_key}_mcce_to_namd.tcl')
            out_lines = self.__write_psfgen_file_for_one_single_case__(case_key)
            with open(out_path, 'w') as fout:
                fout.write(out_lines)
            print(f'The TCL scritp for converting MMCE calculations to NAMD for {case_key} has been deposited at {out_path}.')


    def __write_psfgen_file_for_one_single_case__(self, case_key):
        '''
        Task: Write out the psfgen file needed to build the system for 1 single case
        '''
        #-1st: split protein into chains
        split_pdb_by_chain_lines = f'mol new {self.in_pdb}\n\n'
        pdb_prefix = os.path.basename(self.in_pdb).replace('.pdb', '').strip()
        pdb_dir = os.path.dirname(self.in_pdb)
        for chain in self.chains:
            split_pdb_by_chain_lines += f'set sel [atomselect top "chain {chain}"]\n$sel writepdb {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}.pdb\n\n'

        #-2nd: rename PDB by chain
        rename_lines = ''
        for chain in self.chains:
            rename_lines += f'mol new {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}.pdb\n'
            for line in self.actions_by_chains[case_key][chain]:
                psfgen_priority = line[1]
                if psfgen_priority == 1:
                    rename_lines += line[2] + '\n'
            rename_lines += f'''set sel [atomselect top "all"]
$sel writepdb {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}_renamed.pdb
\n'''

        #-3rd: build psf & patches for each chain
        psf_lines = ''
        for chain in self.chains:
            psf_lines += f'''package require psfgen
resetpsf\n\n'''
            if not self.in_topologies:
                pass
            else:
                for topo in self.in_topologies:
                    psf_lines += f'topology {topo}\n'
            psf_lines += '\n'
            psf_lines += f'segment SER{chain} ' + '{\n' +\
            f'  pdb {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}_renamed.pdb\n'
            for line in self.actions_by_chains[case_key][chain]:
                psfgen_priority = line[1]
                if psfgen_priority == 2 or psfgen_priority == 3:
                    psf_lines += '  ' + line[2] +'\n'
            psf_lines += '}\n\n'
            for line in self.actions_by_chains[case_key][chain]:
                psfgen_priority = line[1]
                if psfgen_priority == 4:
                    psf_lines += '  ' + line[2] +'\n'
            psf_lines += f'''
pdbalias atom ILE CD1 CD
pdbalias atom HOH O OH2
pdbalias residue HOH TIP3
coordpdb {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}_renamed.pdb SER{chain}

guesscoord

writepsf {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}_psfgen.psf
writepdb {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}_psfgen.pdb

'''
        #-4th: merge PSF & PDB from different chains
        merge_lines = ''
        if len(self.chains) == 1:
            merge_lines += f'''mol load psf {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}_psfgen.psf pdb {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}_psfgen.pdb
set al [atomselect top "all"]
$al writepdb {pdb_dir}/{pdb_prefix}_{case_key}_namd_ready.pdb
$al writepsf {pdb_dir}/{pdb_prefix}_{case_key}_namd_ready.psf
'''
        elif len(self.chains) > 1:
            merge_lines += '''package require psfgen
resetpsf

'''
            for chain in self.chains:
                merge_lines += f'readpsf  {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}_psfgen.psf pdb {pdb_dir}/{pdb_prefix}_{case_key}_chain_{chain}_psfgen.pdb\n'
            merge_lines += f'''writepsf {pdb_dir}/{pdb_prefix}_{case_key}_namd_ready.psf
writepdb {pdb_dir}/{pdb_prefix}_{case_key}_namd_ready.pdb'''

        #-5th combine all lines
        psfgen_lines = split_pdb_by_chain_lines + rename_lines + psf_lines + merge_lines
        return psfgen_lines


    def __sort_actions_by_chains__(self, actions_dict: dict=None):
        '''
        Task: Group message for the same chain into 1 group.

        !! May NOT need that as VMD psfgen should be able to build psf for PDB with multiple chians
        as long as each chain is being loaded separately.
        1) just output each chain separately [read Catherine's notes;
        2) then load the corresponding PDB for psfgen during segment creation.
        '''
        if not actions_dict:
            print(f'No action is needed/provided for the current case.')
            return None
        else:
            for case_key in actions_dict.keys():
                lines_by_chains = defaultdict(lambda: [])
                in_lines = actions_dict[case_key]
                for chain in self.chains:
                    for line in in_lines:
                        chain_id = line[0][-1]
                        if chain_id == chain:
                            lines_by_chains[chain].append(line)
                    lines_by_chains[chain] = sorted(lines_by_chains[chain], key=lambda elm: elm[1])
                actions_dict[case_key] = lines_by_chains
        return actions_dict


    def __get_chain_and_residue_types__(self):
        '''
        Task: Get the type of residues involved in MCCE calculations.
        '''
        residues_involved = self.in_ms_df.columns[:-3]
        residue_types_involved = list(set(
            [key[:3] for key in residues_involved]
        ))
        self.residues_involved = residues_involved
        self.residue_types_involved = residue_types_involved
        self.chains = list(set([res[3] for res in self.residues_involved]))


    def __get_per_type_charge_states__(self):
        '''
        Task: Get the charge states involved for each residue type
        present in MCCE calculations.
        '''
        crg_dict = defaultdict(lambda: [])
        for res in self.residues_involved:
            resname = res[:3]
            crg_dict[resname].extend(self.in_ms_df[res].to_list())
        for resname in crg_dict.keys():
            crg_dict[resname] = list(set(crg_dict[resname]))
        self.charge_states_involved = crg_dict


    def __get_message_and_psfgen_lines_for_a_row__(self, in_single_row_df):
        '''
        Task: Generate action messages and TCL lines based on a single input row
        from a larger CSV table containing many MCCE calculations. Each row consists
        of many columns, depending of the number of residues being considered.

        Args:
            in_single_row_df:
                - A single row of self.sel_ms_df being casted in a form of dataFrame,
                following the format of the input example.
        '''
        actions = []
        for res in self.residues_involved:
            resname = res[:3]
            chain_id = res[3]
            resid = int(res[4:-1])
            charge_state = in_single_row_df[res].item()
            segment_id = 'SEG' + str(chain_id)
            segment_id, psfgen_priority, psfgen_line, message_line = \
            self.__write_message_and_psfgen_lines__(
                resname, chain_id, resid, charge_state, segment_id
            )
            actions.append((segment_id, psfgen_priority, psfgen_line, message_line))
        return actions


    def __write_message_and_psfgen_lines__(
            self,
            resname,
            chain_id,
            resid,
            charge_state,
            segment_id,
        ):
        '''
        Task: Generate messages and TCL lines based on info for the input residue.
        The messages are for users to take action if they are to be done manually.
        The lines are for creating a corresponding psfgen TCL file that will build the
        needed NAMD system through VMD.
        '''
        patch_to_use = self.Info.charmm_dict_for_MCCE[resname][str(charge_state)]
        patch_key = patch_to_use[5:].strip()
        if resname == 'CTR':
            message_line = f'Apply patch "{patch_to_use}" to the C-terminal (resid {resid}) of chain {chain_id}.'
            psfgen_line = f'last {patch_key}'
            psfgen_priority = 3
        elif resname == 'NTR':
            #--This function does not work yet because NTR is not yet included in our dict.
            #message_line = f'Apply patch "{patch_to_use}" to the N-terminal (resid {resid}) of chain {chain_id}.'
            #psfgen_line = f'first {patch_key}'
            #message_line = None
            message_line = 'Not implemented.'
            #psfgen_line = None
            psfgen_line = 'Not implemented.'
            psfgen_priority = 2
        else:
            # assuming the resname has been documented as either a regular residue or a patch
            if 'RESI' in patch_to_use:
                #--Regular residue
                message_line = f'Name/Rename residue with resid {resid} of chain {chain_id} to {patch_key}.'
                psfgen_line = f'''set sel [atomselect top "resid {resid} and chain {chain_id}"]
$sel set resname {patch_key}'''
                psfgen_priority = 1
            elif 'PRES' in patch_to_use:
                message_line = f'Apply patch ({patch_to_use}) to resid {resid} of chain {chain_id}.'
                psfgen_line = f'patch {patch_key} {segment_id}:{resid}'
                psfgen_priority = 4
            else:
                #--not implemented
                pass
        return segment_id, psfgen_priority, psfgen_line, message_line
