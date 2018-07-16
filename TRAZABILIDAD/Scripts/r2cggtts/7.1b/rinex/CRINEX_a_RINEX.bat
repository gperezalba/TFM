echo =================================================
echo EJECUTANDO HATANAKA EN FICHEROS RINEX *.yyd
echo =================================================

FOR %%j IN (*.??d) DO (
                        crx2rnx %%j -f
                        echo %%j
                        DEL %%j
                        )