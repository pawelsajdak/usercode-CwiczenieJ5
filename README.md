#Instrukcja poczatkowa cwiczenia

#logowanie na komputer zcoble z fizyk2
ssh -X zcobl.hep

#Ustawianie srodowiska CMSSW dla bash (można dopisc do .bashrc)
export VO_CMS_SW_DIR=/opt/soft/cms
source $VO_CMS_SW_DIR/cmsset_default.sh

#sprawdzenie dostepnych wersji CMSSW
scram list CMSSW

#Inicjalizacja przestrzeni roboczej w wersji CMSSW_8_0_26
cmsrel CMSSW_8_0_26
cd CMSSW_8_0_26/src
cmsenv

# Sciagniecie cwiczenia i kompilacja
git clone https://github.com/konec/usercode-CwiczenieJ5.git UserCode/Cwiczenie
scram b -r

# Puszczenie jobu testowego i zrobienie histogrmau kontrolnego
cd UserCode/Cwiczenie/test
cmsRun cwiczenie.py
./test.py
