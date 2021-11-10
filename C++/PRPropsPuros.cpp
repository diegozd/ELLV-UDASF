#include "PRPropsPuros.hpp"

// Constructor da classe
PRPropsPuros::PRPropsPuros(std::string path) : MGProps(path)
{
    unsigned int inicio = path.rfind("/", path.size()-2) + 1;
    unsigned int tamanho = path.size()-inicio-1;
    std::string nomeCorrente = path.substr(inicio,tamanho);

    CalcPropsPuros();
    
    if (nomeCorrente == "C3" || nomeCorrente == "c3")
    {
        //propano
        mol[0].SOLV = 1;
        mol[0].Tm = 85.46; //Normal melting point (Tm [K])
        mol[0].Tb = 231.11; // Normal boiling point (Tb [K])
        mol[0].Tc = 369.82; //Critical temperature (Tc [K])
        mol[0].Pc = 42.49; //Critical pressure (Pc [bar])
        mol[0].Vc = 202.9;  //Critical Volum (Vc [cm³/mol])
        mol[0].Gf = -24.3;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[0].Hf = -104.7;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[0].Hv = 18.8;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[0].Hfus = 3.524;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[0].VmolRckt =  (mol[0].Pc/(R*mol[0].Tc)*(1/pow(((mol[0].Pc*mol[0].Vc/1000000)/(R*mol[0].Tc)), (1 + pow((1 - Tref/mol[0].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[0].w = 0.152; // fator acêntrico de cada molécula - artigo de 2019
        mol[0].MMw = 44;
    }
    else if (nomeCorrente == "GLP" || nomeCorrente == "Glp" || nomeCorrente == "glp")
    {        
        // 0 - propano
        mol[0].SOLV = 1;
        mol[0].Tm = 85.46; //Normal melting point (Tm [K])
        mol[0].Tb = 231.0; // Normal boiling point (Tb [K])
        mol[0].Tc = 369.9; //Critical temperature (Tc [K])
        mol[0].Pc = 42.5666; //Critical pressure (Pc [bar])
        mol[0].Vc = 200;  //Critical Volum (Vc [cm³/mol])
        mol[0].Gf = -24.3;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[0].Hf = -104.7/1000;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[0].Hv = 18.8;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[0].Hfus = 3.524;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[0].VmolRckt =  (mol[0].Pc/(R*mol[0].Tc)*(1/pow(((mol[0].Pc*mol[0].Vc/1000000)/(R*mol[0].Tc)), (1 + pow((1 - Tref/mol[0].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[0].w = 0.1524; // fator acêntrico de cada molécula - artigo de 2019
        mol[0].MMw = 44.1;

        // 1 - propeno
        mol[1].SOLV = 1;
        mol[1].Tm = 87.9; //Normal melting point (Tm [K])
        mol[1].Tb = 225.43; // Normal boiling point (Tb [K])
        mol[1].Tc = 365;//364.76; //Critical temperature (Tc [K])
        mol[1].Pc = 46.2041;//46.13; //Critical pressure (Pc [bar])
        mol[1].Vc = 181.0;  //Critical Volum (Vc [cm³/mol])
        mol[1].Gf = 62.2;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[1].Hf = 19.7;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[1].Hv = 18.49;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[1].Hfus = 3.003;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[1].VmolRckt =  (mol[1].Pc/(R*mol[1].Tc)*(1/pow(((mol[1].Pc*mol[1].Vc/1000000)/(R*mol[1].Tc)), (1 + pow((1 - Tref/mol[1].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[1].w = 0.1480;//0.142; // fator acêntrico de cada molécula - artigo de 2019
        mol[1].MMw = 42.081;

        // 2 - butano
        mol[2].SOLV = 1;
        mol[2].Tm = 134.86; //Normal melting point (Tm [K])
        mol[2].Tb = 272.6; // Normal boiling point (Tb [K])
        mol[2].Tc = 425.2; //Critical temperature (Tc [K])
        mol[2].Pc = 37.9662; //Critical pressure (Pc [bar])
        mol[2].Vc = 255;  //Critical Volum (Vc [cm³/mol])
        mol[2].Gf = -15.9;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[2].Hf = -126.8;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[2].Hv = 22.44;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[2].Hfus = 4.660;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[2].VmolRckt =  (mol[2].Pc/(R*mol[2].Tc)*(1/pow(((mol[2].Pc*mol[2].Vc/1000000)/(R*mol[2].Tc)), (1 + pow((1 - Tref/mol[2].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[2].w = 0.2010; // fator acêntrico de cada molécula - artigo de 2019
        mol[2].MMw = 58.12;

        // 3 - Buteno
        mol[3].SOLV = 1;
        mol[3].Tm = 87.8; //Normal melting point (Tm [K])
        mol[3].Tb = 266.9; // Normal boiling point (Tb [K])
        mol[3].Tc = 419.60;//419.59; //Critical temperature (Tc [K])
        mol[3].Pc = 40.226;//40.20; //Critical pressure (Pc [bar])
        mol[3].Vc = 239.9;  //Critical Volum (Vc [cm³/mol])
        mol[3].Gf = 70.4;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[3].Hf = -0.5;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[3].Hv = 22.44;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[3].Hfus = 3.849;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[3].VmolRckt =  (mol[3].Pc/(R*mol[3].Tc)*(1/pow(((mol[3].Pc*mol[3].Vc/1000000)/(R*mol[3].Tc)), (1 + pow((1 - Tref/mol[3].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[3].w = 0.187; // fator acêntrico de cada molécula - artigo de 2019
        mol[3].MMw = 56.108;

        // 4 - Cis-buteno-2
        mol[4].SOLV = 1;
        mol[4].Tm = 133.8; //Normal melting point (Tm [K])
        mol[4].Tb = 276.87; // Normal boiling point (Tb [K])
        mol[4].Tc = 435.58; //Critical temperature (Tc [K])
        mol[4].Pc = 42.0577; //Critical pressure (Pc [bar])
        mol[4].Vc = 234.0;  //Critical Volum (Vc [cm³/mol])
        mol[4].Gf = 65.5;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[4].Hf = -7.4;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[4].Hv = 23.43;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[4].Hfus = 7.318;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[4].VmolRckt =  (mol[4].Pc/(R*mol[4].Tc)*(1/pow(((mol[4].Pc*mol[4].Vc/1000000)/(R*mol[4].Tc)), (1 + pow((1 - Tref/mol[4].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[4].w = 0.203; // fator acêntrico de cada molécula - artigo de 2019
        mol[4].MMw = 56.108;

        // 5 - Transbuteno-2
        mol[5].SOLV = 1;
        mol[5].Tm = 167.6; //Normal melting point (Tm [K])
        mol[5].Tb = 274.03; // Normal boiling point (Tb [K])
        mol[5].Tc = 428.63; //Critical temperature (Tc [K])
        mol[5].Pc = 41.0235; //Critical pressure (Pc [bar])
        mol[5].Vc = 238.2;  //Critical Volum (Vc [cm³/mol])
        mol[5].Gf = 63.3;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[5].Hf = -11;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[5].Hv = 22.91;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[5].Hfus = 9.760;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[5].VmolRckt =  (mol[5].Pc/(R*mol[5].Tc)*(1/pow(((mol[5].Pc*mol[5].Vc/1000000)/(R*mol[5].Tc)), (1 + pow((1 - Tref/mol[5].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[5].w = 0.2182; // fator acêntrico de cada molécula - artigo de 2019
        mol[5].MMw = 56.108;

        // 6 - Isobutano
        mol[6].SOLV = 1;
        mol[6].Tm = 113.7; //Normal melting point (Tm [K])
        mol[6].Tb = 261.43; // Normal boiling point (Tb [K])
        mol[6].Tc = 407.98;//408.14; //Critical temperature (Tc [K])
        mol[6].Pc = 36.5491;//36.48; //Critical pressure (Pc [bar])
        mol[6].Vc = 262.7;  //Critical Volum (Vc [cm³/mol])
        mol[6].Gf = -21.4;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[6].Hf = -135;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[6].Hv = 21.4;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[6].Hfus = 4.540;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[6].VmolRckt =  (mol[6].Pc/(R*mol[6].Tc)*(1/pow(((mol[6].Pc*mol[6].Vc/1000000)/(R*mol[6].Tc)), (1 + pow((1 - Tref/mol[6].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[6].w = 0.1812;//0.177; // fator acêntrico de cada molécula - artigo de 2019
        mol[6].MMw = 58.123;

        // 7 - Isobuteno
        mol[7].SOLV = 1;
        mol[7].Tm = 132.4; //Normal melting point (Tm [K])
        mol[7].Tb = 266.25; // Normal boiling point (Tb [K])
        mol[7].Tc = 417.90; //Critical temperature (Tc [K])
        mol[7].Pc = 40.0233;//39.99; //Critical pressure (Pc [bar])
        mol[7].Vc = 238.9;  //Critical Volum (Vc [cm³/mol])
        mol[7].Gf = 58.2;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[7].Hf = -17.1;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[7].Hv = 22.21;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[7].Hfus = 5.920;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[7].VmolRckt =  (mol[7].Pc/(R*mol[7].Tc)*(1/pow(((mol[7].Pc*mol[7].Vc/1000000)/(R*mol[7].Tc)), (1 + pow((1 - Tref/mol[7].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[7].w = 0.1900;//0.189; // fator acêntrico de cada molécula - artigo de 2019
        mol[7].MMw = 56.108;
    }
    else if (nomeCorrente == "C7" || nomeCorrente == "c7")
    {
        //n-heptano
        mol[0].SOLV = 1;
        mol[0].Tm = 182.57; //Normal melting point (Tm [K])
        mol[0].Tb = 371.58; // Normal boiling point (Tb [K])
        mol[0].Tc = 540.26; //Critical temperature (Tc [K])
        mol[0].Pc = 27.36; //Critical pressure (Pc [bar])
        mol[0].Vc = 431.9;  //Critical Volum (Vc [cm³/mol])
        mol[0].Gf = 8.2;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[0].Hf = -657.39;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[0].Hv = 31.73;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[0].Hfus = 14.030;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[0].VmolRckt =  (mol[0].Pc/(R*mol[0].Tc)*(1/pow(((mol[0].Pc*mol[0].Vc/1000000)/(R*mol[0].Tc)), (1 + pow((1 - Tref/mol[0].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[0].w = 0.351; // fator acêntrico de cada molécula - artigo de 2019
        mol[0].MMw = 100.2;
    }
    else if (nomeCorrente == "BZ_mCyC6")
    {
        //benzeno
        mol[0].SOLV = 1;
        mol[0].Tm = 0; //Normal melting point (Tm [K])
        mol[0].Tb = 353.24; // Normal boiling point (Tb [K])
        mol[0].Tc = 562.1; //Critical temperature (Tc [K])
        mol[0].Pc = 49.2439; //Critical pressure (Pc [bar])
        mol[0].Vc = 260;  //Critical Volum (Vc [cm³/mol])
        mol[0].Gf = 0;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[0].Hf = 0;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[0].Hv = 0;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[0].Hfus = 0;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[0].VmolRckt =  (mol[0].Pc/(R*mol[0].Tc)*(1/pow(((mol[0].Pc*mol[0].Vc/1000000)/(R*mol[0].Tc)), (1 + pow((1 - Tref/mol[0].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[0].w = 0.2150; // fator acêntrico de cada molécula - artigo de 2019
        mol[0].MMw = 78.11;

        //methyl-ciclhexane
        mol[1].SOLV = 1;
        mol[1].Tm = 0; //Normal melting point (Tm [K])
        mol[1].Tb = 374.08; // Normal boiling point (Tb [K])
        mol[1].Tc = 572.1; //Critical temperature (Tc [K])
        mol[1].Pc = 34.7537; //Critical pressure (Pc [bar])
        mol[1].Vc = 368;  //Critical Volum (Vc [cm³/mol])
        mol[1].Gf = 0;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[1].Hf = 0;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[1].Hv = 0;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[1].Hfus = 0;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[1].VmolRckt =  (mol[0].Pc/(R*mol[0].Tc)*(1/pow(((mol[0].Pc*mol[0].Vc/1000000)/(R*mol[0].Tc)), (1 + pow((1 - Tref/mol[0].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[1].w = 0.2330; // fator acêntrico de cada molécula - artigo de 2019
        mol[1].MMw = 98.19;
    }
    else if (nomeCorrente == "BZ_C3O")
    {
        //benzeno
        mol[0].SOLV = 1;
        mol[0].Tm = 0; //Normal melting point (Tm [K])
        mol[0].Tb = 353.24; // Normal boiling point (Tb [K])
        mol[0].Tc = 562.1; //Critical temperature (Tc [K])
        mol[0].Pc = 49.2439; //Critical pressure (Pc [bar])
        mol[0].Vc = 260;  //Critical Volum (Vc [cm³/mol])
        mol[0].Gf = 0;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[0].Hf = 0;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[0].Hv = 0;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[0].Hfus = 0;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[0].VmolRckt =  (mol[0].Pc/(R*mol[0].Tc)*(1/pow(((mol[0].Pc*mol[0].Vc/1000000)/(R*mol[0].Tc)), (1 + pow((1 - Tref/mol[0].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[0].w = 0.2150; // fator acêntrico de cada molécula - artigo de 2019
        mol[0].MMw = 78.11;

        //2-propanol
        mol[1].SOLV = 1;
        mol[1].Tm = 0; //Normal melting point (Tm [K])
        mol[1].Tb = 386.75; // Normal boiling point (Tb [K])
        mol[1].Tc = 570; //Critical temperature (Tc [K])
        mol[1].Pc = 65.2; //Critical pressure (Pc [bar])
        mol[1].Vc = 176;  //Critical Volum (Vc [cm³/mol])
        mol[1].Gf = 0;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[1].Hf = 0;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[1].Hv = 0;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[1].Hfus = 0;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[1].VmolRckt =  (mol[0].Pc/(R*mol[0].Tc)*(1/pow(((mol[0].Pc*mol[0].Vc/1000000)/(R*mol[0].Tc)), (1 + pow((1 - Tref/mol[0].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[1].w = 0.5546; // fator acêntrico de cada molécula - artigo de 2019
        mol[1].MMw = 56.06;
    }
    else if (nomeCorrente == "BZ_mCyC6_H2O")
    {
        //benzeno
        mol[0].SOLV = 1;
        mol[0].Tm = 0; //Normal melting point (Tm [K])
        mol[0].Tb = 353.24; // Normal boiling point (Tb [K])
        mol[0].Tc = 562.1; //Critical temperature (Tc [K])
        mol[0].Pc = 49.2439; //Critical pressure (Pc [bar])
        mol[0].Vc = 260;  //Critical Volum (Vc [cm³/mol])
        mol[0].Gf = 0;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[0].Hf = 0;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[0].Hv = 0;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[0].Hfus = 0;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[0].VmolRckt =  (mol[0].Pc/(R*mol[0].Tc)*(1/pow(((mol[0].Pc*mol[0].Vc/1000000)/(R*mol[0].Tc)), (1 + pow((1 - Tref/mol[0].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[0].w = 0.2150; // fator acêntrico de cada molécula - artigo de 2019
        mol[0].MMw = 78.11;

        //methyl-ciclhexane
        mol[1].SOLV = 1;
        mol[1].Tm = 0; //Normal melting point (Tm [K])
        mol[1].Tb = 374.08; // Normal boiling point (Tb [K])
        mol[1].Tc = 572.1; //Critical temperature (Tc [K])
        mol[1].Pc = 34.7537; //Critical pressure (Pc [bar])
        mol[1].Vc = 368;  //Critical Volum (Vc [cm³/mol])
        mol[1].Gf = 0;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[1].Hf = 0;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[1].Hv = 0;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[1].Hfus = 0;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[1].VmolRckt =  (mol[0].Pc/(R*mol[0].Tc)*(1/pow(((mol[0].Pc*mol[0].Vc/1000000)/(R*mol[0].Tc)), (1 + pow((1 - Tref/mol[0].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[1].w = 0.2330; // fator acêntrico de cada molécula - artigo de 2019
        mol[1].MMw = 98.19;

        //água
        mol[2].SOLV = 1;
        mol[2].Tm = 0; //Normal melting point (Tm [K])
        mol[2].Tb = 373.15; // Normal boiling point (Tb [K])
        mol[2].Tc = 647.3; //Critical temperature (Tc [K])
        mol[2].Pc = 221.2; //Critical pressure (Pc [bar])
        mol[2].Vc = 57.1;  //Critical Volum (Vc [cm³/mol])
        mol[2].Gf = 0;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[2].Hf = 0;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[2].Hv = 0;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[2].Hfus = 0;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[2].VmolRckt =  (mol[0].Pc/(R*mol[0].Tc)*(1/pow(((mol[0].Pc*mol[0].Vc/1000000)/(R*mol[0].Tc)), (1 + pow((1 - Tref/mol[0].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[2].w = 0.3440; // fator acêntrico de cada molécula - artigo de 2019
        mol[2].MMw = 18.02;
    }
    else if (nomeCorrente == "BZ_nC6_H2O")
    {
        //benzeno
        mol[0].SOLV = 1;
        mol[0].Tm = 0; //Normal melting point (Tm [K])
        mol[0].Tb = 353.24; // Normal boiling point (Tb [K])
        mol[0].Tc = 562.1; //Critical temperature (Tc [K])
        mol[0].Pc = 49.2439; //Critical pressure (Pc [bar])
        mol[0].Vc = 260;  //Critical Volum (Vc [cm³/mol])
        mol[0].Gf = 0;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[0].Hf = 0;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[0].Hv = 0;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[0].Hfus = 0;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[0].VmolRckt =  (mol[0].Pc/(R*mol[0].Tc)*(1/pow(((mol[0].Pc*mol[0].Vc/1000000)/(R*mol[0].Tc)), (1 + pow((1 - Tref/mol[0].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[0].w = 0.2150; // fator acêntrico de cada molécula - artigo de 2019
        mol[0].MMw = 78.11;

        //n-hexano
        mol[1].SOLV = 1;
        mol[1].Tm = 0; //Normal melting point (Tm [K])
        mol[1].Tb = 341.88; // Normal boiling point (Tb [K])
        mol[1].Tc = 507.9; //Critical temperature (Tc [K])
        mol[1].Pc = 29.8807; //Critical pressure (Pc [bar])
        mol[1].Vc = 368;  //Critical Volum (Vc [cm³/mol])
        mol[1].Gf = 0;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[1].Hf = 0;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[1].Hv = 0;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[1].Hfus = 0;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[1].VmolRckt =  (mol[0].Pc/(R*mol[0].Tc)*(1/pow(((mol[0].Pc*mol[0].Vc/1000000)/(R*mol[0].Tc)), (1 + pow((1 - Tref/mol[0].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[1].w = 0.3007; // fator acêntrico de cada molécula - artigo de 2019
        mol[1].MMw = 86.18;

        //água
        mol[2].SOLV = 1;
        mol[2].Tm = 0; //Normal melting point (Tm [K])
        mol[2].Tb = 373.15; // Normal boiling point (Tb [K])
        mol[2].Tc = 647.3; //Critical temperature (Tc [K])
        mol[2].Pc = 221.2; //Critical pressure (Pc [bar])
        mol[2].Vc = 57.1;  //Critical Volum (Vc [cm³/mol])
        mol[2].Gf = 0;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[2].Hf = 0;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[2].Hv = 0;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[2].Hfus = 0;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[2].VmolRckt =  (mol[0].Pc/(R*mol[0].Tc)*(1/pow(((mol[0].Pc*mol[0].Vc/1000000)/(R*mol[0].Tc)), (1 + pow((1 - Tref/mol[0].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[2].w = 0.3440; // fator acêntrico de cada molécula - artigo de 2019
        mol[2].MMw = 18.02;
    }
    else if (nomeCorrente == "BZ_nC6")
    {
        //benzeno
        mol[0].SOLV = 1;
        mol[0].Tm = 0; //Normal melting point (Tm [K])
        mol[0].Tb = 353.24; // Normal boiling point (Tb [K])
        mol[0].Tc = 562.1; //Critical temperature (Tc [K])
        mol[0].Pc = 49.2439; //Critical pressure (Pc [bar])
        mol[0].Vc = 260;  //Critical Volum (Vc [cm³/mol])
        mol[0].Gf = 0;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[0].Hf = 0;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[0].Hv = 0;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[0].Hfus = 0;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[0].VmolRckt =  (mol[0].Pc/(R*mol[0].Tc)*(1/pow(((mol[0].Pc*mol[0].Vc/1000000)/(R*mol[0].Tc)), (1 + pow((1 - Tref/mol[0].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[0].w = 0.2150; // fator acêntrico de cada molécula - artigo de 2019
        mol[0].MMw = 78.11;

        //n-hexano
        mol[1].SOLV = 1;
        mol[1].Tm = 0; //Normal melting point (Tm [K])
        mol[1].Tb = 341.88; // Normal boiling point (Tb [K])
        mol[1].Tc = 507.9; //Critical temperature (Tc [K])
        mol[1].Pc = 29.8807; //Critical pressure (Pc [bar])
        mol[1].Vc = 368;  //Critical Volum (Vc [cm³/mol])
        mol[1].Gf = 0;  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
        mol[1].Hf = 0;  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
        mol[1].Hv = 0;  //Standard enthalpy of vaporization at Tb K (Hv [kJ/mol])
        mol[1].Hfus = 0;  //Standard enthalpy of fusion (Hfus [kJ/mol])
        mol[1].VmolRckt =  (mol[0].Pc/(R*mol[0].Tc)*(1/pow(((mol[0].Pc*mol[0].Vc/1000000)/(R*mol[0].Tc)), (1 + pow((1 - Tref/mol[0].Tc), 2/7)))))/1000; // Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
        mol[1].w = 0.3007; // fator acêntrico de cada molécula - artigo de 2019
        mol[1].MMw = 86.18;
    }


    CalcTbPR(1.01325);
}


void PRPropsPuros::CalcTbPR(double P)    // P em bar
{
    double T, FT, FTdT, FlT;
    double dT = 0.01;
    double e = 0.0001;
    unsigned int iter = 0;
    TbPR.resize(nMol);
    for (unsigned int i = 0 ; i < TbPR.size() ; i++)
    {
        T = mol[i].Tb;
        iter = 0;
        do
        {
            FT = FTPR(T,P,i);
            FTdT = FTPR(T+dT,P,i);
            FlT = (FTdT - FT)/dT;
            T = T - FT/FlT;
            iter++;
            // if (iter >= maxiter) std::cout << "CalcTbPR atingiu maxiter na molecula " << i << std::endl;
        } while ((abs(FT) > e) && (iter < maxiter));
        if (isnan(T)) T = mol[i].Tb;
        TbPR[i] = T;
    }
}

double PRPropsPuros::FZPR(double Z, double A, double B)
{
    double FZ, FZdZ, FlZ;
    double dZ = 0.0001;
    double e = 0.00001;
    unsigned int iter = 0;
    do
    {
        FZ = pow(Z,3) - (1-B)*pow(Z,2) + (A - 2*B - 3*pow(B,2))*Z - (A*B - pow(B,2) - pow(B,3));
        FZdZ =  pow((Z+dZ),3) - (1-B)*pow((Z+dZ),2) + (A - 2*B - 3*pow(B,2))*(Z+dZ) - (A*B - pow(B,2) - pow(B,3));
        FlZ = (FZdZ - FZ)/dZ;
        Z = Z - FZ/FlZ;
        iter++;
        if (iter >= maxiter) std::cout << "FZPR atingiu maxiter " << std::endl;
    } while ((abs(FZ) > e) && (iter < maxiter));
    if (Z < B) Z = B;
    return Z;
}

double PRPropsPuros::FZLPR(double A, double B, double q)
{
    double Z0 = B;
    double Z;
    double sig = 1+pow(2,0.5);
    double eps = 1-pow(2,0.5);
    double e = 0.00001;
    double diff = 0;
    unsigned int iter = 0;
    do
    {
        Z = B + (Z0 + eps*B)*(Z0 + sig*B)*((1 + B - Z0)/(q*B));
        diff = abs(Z-Z0)/Z;
        Z0 = Z;
        iter++;
        if (iter >= maxiter) std::cout << "FZLPR atingiu maxiter " << std::endl;
    } while ((diff > e) && (iter < maxiter));
    if (Z < B) Z = B;

    return Z;
}

double PRPropsPuros::FZVPR(double A, double B, double q)
{
    double Z0 = 1;
    double Z;
    double sig = 1+pow(2,0.5);
    double eps = 1-pow(2,0.5);
    double e = 0.00001;
    double diff = 0;
    unsigned int iter = 0;
    do
    {
        Z = 1 + B - q*B*((Z0 - B)/((Z0 + eps*B)*(Z0+sig*B)));
        diff = abs(Z-Z0)/Z;
        Z0 = Z;
        iter++;
        // if (iter >= maxiter) std::cout << "FZVPR atingiu maxiter " << Z  << " " << diff << " " << A << " " << B << " " << q << std::endl;
    } while ((diff > e) && (iter < maxiter));
    if (Z < B) Z = B;
    return Z;
}



double PRPropsPuros::FTPR(double T, double P, unsigned int i)
{
    double Zv, Zl, Vl, Vv;
    CtsPRStr K = calcCtesPR(T, P, i);
    
    // Zv = FZPR(1, K.A, K.B);
    // Zl = FZPR(K.B, K.A, K.B);
    Zv = FZVPR(K.A, K.B, K.q);     //vapor
    Zl = FZLPR(K.A, K.B, K.q);            //liquido

    Vv = Zv*R*T/P;
    Vl = Zl*R*T/P;

    return P*(Vv - Vl) - R*T*log((Vv-K.b)/(Vl-K.b)) + (K.a*K.alpha/(2*K.b*pow(2,0.5)))*(log((Vv + K.b*(1-pow(2,0.5)))/(Vv + K.b*(1+pow(2,0.5)))) - log((Vl + K.b*(1-pow(2,0.5)))/(Vl + K.b*(1+pow(2,0.5)))));
}

PRPropsPuros::CtsPRStr PRPropsPuros::calcCtesPR(double T, double P,  unsigned int i)
{
    CtsPRStr out;

    out.Tc = mol[i].Tc;
    out.Pc = mol[i].Pc;
    out.w = mol[i].w;

    out.a = 0.457235529*pow(R*out.Tc,2)/out.Pc;
    out.b = 0.0777960739*R*out.Tc/out.Pc;
    out.Tr = T/out.Tc;
    

    if (out.w <= 0.49) out.kappa = 0.37464 + 1.54226*out.w - 0.269926*pow(out.w,2);
    else out.kappa = 0.379642 + (1.48503 - (0.164423 - 0.016666*out.w)*out.w)*out.w;

    out.alpha = pow((1 + out.kappa*(1-pow(out.Tr,0.5))),2);
    out.a_alpha = out.a * out.alpha;
    out.A = out.a_alpha*P/(pow(R*T,2));
    out.B = out.b*P/(R*T);
    out.q = out.a_alpha/(out.b*R*T);

    // if (isnan(out.a_alpha))  std::cout << " a_alpha:" << out.a_alpha << " a:" << out.a << " alpha:" << out.alpha << " Tr:" << out.Tr << " Tc:" << out.Tc << " T:" << T << std::endl;

    return out;
}

double PRPropsPuros::HRPR(double T, double P, unsigned int i, bool EstFis)
{
    double Z, dadT;
    CtsPRStr K = calcCtesPR(T, P, i);

    // if (EstFis == 1) Z = FZPR(1, K.A, K.B);     //vapor
    // else Z = FZPR(K.B, K.A, K.B);            //liquido
    if (EstFis == 1) Z = FZVPR(K.A, K.B, K.q);     //vapor
    else Z = FZLPR(K.A, K.B, K.q);            //liquido
    
    dadT = -0.457235529*pow(R*K.Tc,2)*K.kappa*pow((K.alpha/(K.Tc*T)),0.5)/K.Pc;

    return Z-1 + ((T*dadT-K.a*K.alpha)/(R*T*2*K.b*pow(2,0.5)))*log((Z+(1+pow(2,0.5))*K.B)/(Z+(1-pow(2,0.5))*K.B));
}

double PRPropsPuros::SRPR(double T, double P, unsigned int i, bool EstFis)
{
    double Z, dadT;
    CtsPRStr K = calcCtesPR(T, P, i);

    // if (EstFis == 1) Z = FZPR(1, K.A, K.B);     //vapor
    // else Z = FZPR(K.B, K.A, K.B);            //liquido
    if (EstFis == 1) Z = FZVPR(K.A, K.B, K.q);     //vapor
    else Z = FZLPR(K.A, K.B, K.q);            //liquido
    
    dadT = -0.457235529*pow(R*K.Tc,2)*K.kappa*pow((K.alpha/(K.Tc*T)),0.5)/K.Pc;

    return log(Z-K.B) + (dadT/(R*2*K.b*pow(2,0.5)))*log((Z+(1+pow(2,0.5))*K.B)/(Z+(1-pow(2,0.5))*K.B));
}

double PRPropsPuros::GRPR(double T, double P, unsigned int i, bool EstFis)
{
    double HR_RT = HRPR(T, P, i, EstFis);
    double SR_R = SRPR(T, P, i, EstFis);
    double GR_RT = HR_RT - SR_R;
    return GR_RT;
}

double PRPropsPuros::FTPR_viaG(double T, double P, unsigned int i)
{
    double GRl, GRv;
    GRv = GRPR(T, P, i, 1);
    GRl = GRPR(T, P, i, 0);
    return GRv - GRl;
}

void PRPropsPuros::CalcTbPR_viaG(double P)    // P em bar
{
    double T, FT, FTdT, FlT;
    double dT = 0.01;
    double e = 0.0001;
    unsigned int iter = 0;

    TbPR_viaG.resize(nMol);
    for (unsigned int i = 0 ; i < TbPR.size() ; i++)
    {
        T = mol[i].Tb;
        iter = 0;
        do
        {
            FT = FTPR_viaG(T,P,i);
            FTdT = FTPR_viaG(T+dT,P,i);
            FlT = (FTdT - FT)/dT;
            T = T - FT/FlT;
            iter++;
            if (iter >= maxiter) std::cout << "CalcTbPR_viaG atingiu maxiter " << std::endl;
        } while ((abs(FT) > e) && (iter < maxiter));
        TbPR_viaG[i] = T;
    }
}

double PRPropsPuros::CalcRho(double T, double P, unsigned int i)    // T em K e P em bar
{
    double Z;
    CtsPRStr K = calcCtesPR(T, P, i);
    // Z = FZPR(K.B, K.A, K.B);
    Z = FZLPR(K.A, K.B, K.q);            //liquido
    return  P*mol[i].MMw/(Z*R*T)/1000;
}