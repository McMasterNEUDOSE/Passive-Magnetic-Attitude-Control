function[ggname,attname,Bname,vname,tname]=fileselec(days)
%base='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/.mat';
if days==1
    ggname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_gravgradbody_1s_21Jun2022_22Jun2022.mat';
    attname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/att_ecf_1s_21Jun2022_22Jun2022.mat';
    Bname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_bodymagdata_1s_21Jun2022_22Jun2022.mat';
    vname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_bodyvel_1s_21Jun2022_22Jun2022.mat';
    tname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/time_1s_21Jun2022_22Jun2022.mat';
elseif days==2
    ggname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_gravgradbody_1s_21Jun2022_23Jun2022.mat';
    attname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/att_ecf_1s_21Jun2022_23Jun2022.mat';
    Bname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_bodymagdata_1s_21Jun2022_23Jun2022.mat';
    vname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_bodyvel_1s_21Jun2022_23Jun2022.mat';
    tname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/time_1s_21Jun2022_23Jun2022.mat'; 
    
elseif days==4
    ggname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_gravgradbody_1s_21Jun2022_25Jun2022.mat';
    attname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/att_ecf_1s_21Jun2022_25Jun2022.mat';
    Bname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_bodymagdata_1s_21Jun2022_25Jun2022.mat';
    vname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_bodyvel_1s_21Jun2022_25Jun2022.mat';
    tname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/time_1s_21Jun2022_25Jun2022.mat';
    
 elseif days==5
    ggname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_gravgradbody_1s_21Jun2022_26Jun2022.mat';
    attname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/att_ecf_1s_21Jun2022_26Jun2022.mat';
    Bname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_bodymagdata_1s_21Jun2022_26Jun2022.mat';
    vname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_bodyvel_1s_21Jun2022_26Jun2022.mat';
    tname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/time_1s_21Jun2022_26Jun2022.mat';

elseif days==7
    ggname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_gravgradbody_1s_21Jun2022_28Jun2022.mat';
    attname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/att_ecf_1s_21Jun2022_28Jun2022.mat';
    Bname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_bodymagdata_1s_21Jun2022_28Jun2022.mat';
    vname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_bodyvel_1s_21Jun2022_28Jun2022.mat';
    tname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/time_1s_21Jun2022_28Jun2022.mat';
 elseif days==10
    ggname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_gravgradbody_1s_21Jun2022_1Jul2022.mat';
    attname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/att_ecf_1s_21Jun2022_1Jul2022.mat';
    Bname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_bodymagdata_1s_21Jun2022_1Jul2022.mat';
    vname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_bodyvel_1s_21Jun2022_1Jul2022.mat';
    tname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/time_1s_21Jun2022_1Jul2022.mat';   
 elseif days==14
    ggname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_gravgradbody_1s_21Jun2022_5Jul2022.mat';
    attname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/att_ecf_1s_21Jun2022_5Jul2022.mat';
    Bname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_bodymagdata_1s_21Jun2022_5Jul2022.mat';
    vname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/ecf_bodyvel_1s_21Jun2022_5Jul2022.mat';
    tname='/Users/ZanBarovier/Documents/McMaster/Neudose/PMACDATA/time_1s_21Jun2022_5Jul2022.mat';
end
end
