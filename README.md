#Instrukcja poczatkowa - wstęp do trygera CMS

#Cwiczenie wykonywane jest ma naszynie wirtualnej dzialajacej na 
#komputerze zcobl w domenie hep. W celu zalozenia konta trzeba sie skontaktowac z A. Kalinowskim lub poprosić o to M. Koneckiego.
#Mając konto usera <useid> oraz password można zalogować się na komputer zcobl.
#logowanie usera <userid> na komputer zcobl z sieci FUW (np. tempac):
ssh -X <userid>@zcobl.hep -p 1921

#Sprawdzić czy zmienna środowiskowa VO_CMS_SW_DIR jest ustawiona, ewentalnie ustawić jak poniżej
#Ustawianie srodowiska CMSSW dla bash (można dopisc do .bashrc):
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

#sprawdzenie dostepnych wersji CMSSW
scram list CMSSW

#Inicjalizacja przestrzeni roboczej w wersji CMSSW_10_6_22
cmsrel CMSSW_10_6_22
cd CMSSW_10_6_22/src
cmsenv

# Sciagniecie cwiczenia i kompilacja
git clone https://github.com/konec/usercode-CwiczenieJ5.git UserCode/Cwiczenie
scram b -r

# Puszczenie jobu testowego i zrobienie histogrmau kontrolnego
cd UserCode/Cwiczenie/test
cmsRun cwiczenie.py
./test.py
