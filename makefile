all: build makevenv	

makevenv:
	python3 -m venv venv
	venv/bin/pip install wheel
	venv/bin/pip install -r PYTHON_REQUIREMENTS.txt
	
build:
	nuget install -PreRelease YC.QuickGraph -OutputDirectory pseudoassembly/packages
	msbuild -property:Configuration=Release alignment/exact_match/Aligner3/Aligner3.csproj
	chmod +x alignment/exact_match/Aligner3/bin/Release/Aligner3.exe
	msbuild -property:Configuration=Release pseudoassembly/SamplerEulerianEfficient/SamplerEulerianEfficient.csproj
	chmod +x pseudoassembly/SamplerEulerianEfficient/bin/Release/SamplerEulerianEfficient.exe
