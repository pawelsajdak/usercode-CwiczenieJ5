#Instrukcja poczatkowa cwiczenia

#ćwiczenie wykonywane jest ma naszynie wirtualnej dzialajacej na 
#komputerze zcobl w domenie hep. W celu zalozenia konta trzeba sie skontaktowac z A. Kalinowskim
#logowanie usera <userid> na komputer zcobl z sieci FUW:
ssh -X <userid>@zcobl.hep -p 2222

#Ustawianie srodowiska CMSSW dla bash (można dopisc do .bashrc)
export VO_CMS_SW_DIR=/opt/soft/cms
source $VO_CMS_SW_DIR/cmsset_default.sh

#sprawdzenie dostepnych wersji CMSSW
scram list CMSSW

#Inicjalizacja przestrzeni roboczej w wersji CMSSW_10_6_4
cmsrel CMSSW_10_6_4
cd CMSSW_10_6_24/src
cmsenv

# Sciagniecie cwiczenia i kompilacja
git clone https://github.com/konec/usercode-CwiczenieJ5.git UserCode/Cwiczenie
scram b -r

# Puszczenie jobu testowego i zrobienie histogrmau kontrolnego
cd UserCode/Cwiczenie/test
cmsRun cwiczenie.py
./test.py
