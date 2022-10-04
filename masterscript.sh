#!/bin/bash
myArray=(1a_import_gtex_hypo.R 1b_map_gtex_hypo.R 1c_hypo_inst.R 2a_vinuelacolocislets.R 2b_compirsoneQTLgtex.R 2c_hyprcolocmap.R 2d_VEP.R 2d_VEP.sh 3a_fibmi_inst_multibi.R 3b_fibmi_harmonisemultibi.R 4_plotfigures.R 5_writetables.R 6_writetext.R)
for str in ${myArray[@]}; do
chmod u+x ./$str
done
echo 'Initializing 1a_import_gtex_hypo.R' && ./1a_import_gtex_hypo.R && echo 'Initializing 1b_map_gtex_hypo.R' && ./1b_map_gtex_hypo.R && echo 'Initializing 1c_hypo_inst.R' && ./1c_hypo_inst.R && echo 'Initializing 2a_vinuelacolocislets.R' && ./2a_vinuelacolocislets.R && echo 'Initializing 2b_compirsoneQTLgtex.R' && ./2b_compirsoneQTLgtex.R && echo 'Initializing 2c_hyprcolocmap.R' && ./2c_hyprcolocmap.R && echo 'Initializing 2d_VEP.R' && ./2d_VEP.R && echo 'Initializing 2d_VEP.sh' && ./2d_VEP.sh && echo 'Initializing 3a_fibmi_inst_multibi.R' && ./3a_fibmi_inst_multibi.R && echo 'Initializing 3b_fibmi_harmonisemultibi.R' && ./3b_fibmi_harmonisemultibi.R && echo 'Initializing 4_plotfigures.R' && ./4_plotfigures.R && echo 'Initializing 5_writetables.R' && ./5_writetables.R && echo 'Initializing 6_writetext.R' && ./6_writetext.R && echo 'The master script finished without errors'
