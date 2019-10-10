## How to build on window (Tan's note)

-	Install Cmake
-	Install Microsoft Visual Studio 2015
-	checkout branch window
-	Install libdivsufsort and zlib using the script `build_win.sh` 
-	Open Hera0.1.1/Hera0.1.1.sln by MS Studio
-	Install pthread using Nudget (https://docs.microsoft.com/en-us/nuget/consume-packages/install-use-packages-visual-studio)
-	Add ../local/lib to additional lib directories (https://stackoverflow.com/questions/4445418/how-to-add-additional-libraries-to-visual-studio-project)
-	Add ../local/include to addtional include directories (https://support.pixelink.com/support/solutions/articles/3000044961-configuring-visual-studio-for-c-c-projects)
-	Ctrl + Shift + B to build (remember to build Release instead of Debug)
-	target file will be available at `Hera0.1.1/x64/Release`
-	copy the dll files from local/bin into the Release folder above
-	enjoy it

