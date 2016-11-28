# nTupleMaker

cmsrel CMSSW_8_0_24
cd src
mkdir UserCode
cd Usercode
git clone git@github.com:halilg/nTupleMaker.git
cd nTupleMaker
cmsenv
scram b
./do_test.zsh
