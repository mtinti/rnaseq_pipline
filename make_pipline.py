# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 11:52:11 2014

@author: mtinti-x
"""
import os

replace_list = eval(open('vars.txt').read())
template_R = open('_template_count.R').read()
template_mqc = open('_template_align_barcode.sh').read()

#print(template_file)
if __name__ == '__main__':
    run_all_content = ''
    for dictionary in replace_list:

        sh_script_name = 'j'+dictionary['base_fastq']+'.sh'
        run_all_content+='chmod +x '+sh_script_name+'\n'
        run_all_content+='qsub '+sh_script_name+'\n'
        run_all_content+='mv '+sh_script_name+' '+dictionary['experiment']+'\n'
        
        template_file = open('_template_align_barcode.sh').read()
        sh_script_content = template_file.format(
            g_version=dictionary['g_version'],
            base_fastq=dictionary['base_fastq'],
            experiment=dictionary['experiment'],
            )

        path_to_bam = os.path.join(dictionary['experiment'],'data',dictionary['base_fastq'])
        file_list = ['.sorted.bam']
        
        r_count = template_R.format(
        count_file = os.path.join(dictionary['experiment'], dictionary['base_fastq']+'count.txt'),
        gtf_file = os.path.join('genomes', dictionary['g_version'], dictionary['g_version']+'.gtf'),
        bam_files = ',\n'.join([ '\"'+os.path.join(path_to_bam, dictionary['base_fastq']+n)+'\"' for n in file_list ])
        )

        open('count_'+dictionary['base_fastq']+'.R','w').write(r_count)
        open(sh_script_name,'w').write(sh_script_content)
        
        
  

    run_all_content+='mv '+'run_all_'+dictionary['experiment']+'.sh'+' '+dictionary['experiment']+'\n'
    open('run_all_'+dictionary['experiment']+'.sh','w').write(run_all_content)
    
    print(replace_list)
    template_R = open('_template_count.R').read()
    r_count_all = template_R.format(
        count_file = os.path.join(dictionary['experiment'],dictionary['experiment']+'_count_all.txt'),
        gtf_file = os.path.join('genomes', dictionary['g_version'], dictionary['g_version']+'.gtf'),
        bam_files = ',\n'.join([ '\"'+os.path.join(dictionary['experiment'],'data',dictionary['base_fastq'],dictionary['base_fastq']+'.sorted.bam')+'\"' for dictionary in replace_list ])
        )
    open('count_all_'+dictionary['experiment']+'.R','w').write(r_count_all)
    
    template_mqc = open('_template_multiqc.sh').read()
    mqc = template_mqc.format(experiment=replace_list[0]['experiment'])
    open('mqc_all_'+replace_list[0]['experiment']+'.sh','w').write(mqc)
