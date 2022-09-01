# 101_Dalmatians
Brief description of the directories:<br>
<br>
##########################################<br>
<i>Computational_Models/Full_Models/Original/</i><br>
Sets of features extracted through computational modeling from the movie. Movie-related features fall into two categories: 
i) low-level acoustic and visual features; and ii) high-level semantic descriptors. 
The first category comprises fine-grained features extracted from the auditory (e.g., spectral and sound envelope properties to account for frequency- and amplitude-based modulations) and visual streams (e.g., set of static Gabor-like filters and motion energy information based on their spatiotemporal integration). 
As concerns the set of high-level features, we performed manual annotation of both visual (e.g., a close-up of a dog) and sound-based (e.g., barking of a dog) natural and artificial categories, and we exploited machine learning techniques (e.g., word embedding) to define a group of descriptors representing the semantic properties of the auditory and visual streams. <br>
<br>
##########################################<br>
<i>Computational_Models/Full_Models/Cleaned_from_editing</i><br>
From the above-mentioned computational models we regressed out information related to editing features (e.g., scene transitions, cuts, dialogues, music and audio descriptions).<br>
<br>
##########################################<br>
<i>Computational_Models/Reduced_Models/</i><br>
Selected subsets of acoustic, visual and semantic features, with equal dimensionality and high predictive power.<br>
<br>
##########################################<br>
<i>Results/</i><br>
The directory contains brain masks, ISC onegroup, ISC between groups, Mediation, as well as TRW results.<br>
<br>
##########################################<br>
<i>Code/</i><br>
Matlab functions needed to perform all the analyses<br>
<br>
<br>
More details are available in the preprint: https://www.biorxiv.org/content/10.1101/2022.03.14.484231v1<br>
For fMRI data, please refer to OSF: https://osf.io/j8x6h/<br>

