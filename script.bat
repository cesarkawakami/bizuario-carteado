@FOR %%A in (%*) DO @(
md %%A
echo n=%%A> %%A\Makefile
type Makefile >> %%A\Makefile
rem copy modelo.cpp %%A\%%A.cpp >lixo.txt
)
@del lixo.txt
