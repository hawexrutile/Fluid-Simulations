squeue | grep alan | grep -Eo '[0-9]{6}' | xargs -I {} scancel {}
rm M7-*.cpp