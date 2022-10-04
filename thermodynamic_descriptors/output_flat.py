# -*- coding: utf-8 -*-#

# -------------------------------------------------------------------------------
# Name:         output_flat.py
# Description:
# Author:       Administrator
# Date:         2022/10/4
# -------------------------------------------------------------------------------
from glob import glob
import sys
from ase.calculators.calculator import InputError
from catmap.model import ReactionModel


mkm_job_output = ['production_rate','rate_control','coverage','rate', 'selectivity_control', 'rxn_order']

logfile = glob('*.log')
if len(logfile) > 1:
    raise InputError('Ambiguous logfile. Ensure that only one file ends with .log')
model = ReactionModel(setup_file=logfile[0])

# output_variable = 'production_rate'
for output_variable in mkm_job_output:
    # if output_variable == 'rate_control':
    if output_variable in ['rate_control', "selectivity_control","rxn_order"]:
        dim = 2
    else:
        dim = 1

    labels = model.output_labels[output_variable]

    def flatten_2d(output):
        "Helper function for flattening rate_control output"
        flat = []
        for x in output:
            flat+= x
        return flat

    #flatten rate_control labels
    # if output_variable == 'rate_control':
    if output_variable in ['rate_control',"selectivity_control","rxn_order"]:
        flat_labels = []
        for i in labels[0]:
            for j in labels[1]:
                flat_labels.append('d'+i+'/d'+j)
        labels = flat_labels

    #flatten elementary-step specific labels
    if output_variable in ['rate','rate_constant','forward_rate_constant','reverse_rate_constant']:
        str_labels = []
        for label in labels:
            states = ['+'.join(s) for s in label]
            if len(states) == 2:
                new_label = '<->'.join(states)
            else:
                new_label = states[0]+'<->'+states[1]+'->'+states[2]
            str_labels.append(new_label)
        labels = str_labels
    # print("################### output_variable ",output_variable )
    # print("################### labels ",labels )
    # print("################### model.descriptor_names ",list(model.descriptor_names))
    table = '\t'.join(list(['descriptor-'+d for d in model.descriptor_names])+list(labels))+'\n'

    for pt, output in getattr(model,output_variable+'_map'):
        if dim == 2:
            output = flatten_2d(output)
        table += '\t'.join([str(float(i)) for i in pt+output])+'\n'

    f = open(output_variable+'_table.txt','w')
    f.write(table)
    f.close()