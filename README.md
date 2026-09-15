## Preprocessing:
Download all 4 preprocessing files:  EEGproject_preprocessing.m , avref.eog.mat , channelselection.mat , trialdef.mat  
Open the EEGproject_preprocessing script in Matlab (I used version Matlab R2024b).  
Enter the path to your project folder and your spm path into the EEGproject_preprocessing script:  
```
% Set path
project_root = ['C:\XXX'];
spm_path = ['C:\XXX'];
```  
Make sure your folder structure looks EXACTLY like this:
```
    data /
        00Behavioural /
            images
            logs
            neuronavigation /
                Gian_ID04.sfp
        01EEG /
            raw /
                SPNCartoons_ID04.bdf
            spm /
                channelselection.mat
                avref_eog.mat
                trialdef.mat
        05Anat
```
Then you can run the script.
