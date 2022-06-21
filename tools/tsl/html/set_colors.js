function setColors(radio)  { 
    var scheme = radio.value;
    var name;

    for(var i=0; i<15; i++)  {
        name = 'symbols'+i
        document.form.elements[name].value="";
        name = 'color'+i
        document.form.elements[name].value="";
    }



    if (scheme=="amino_bw")  {
                document.form.symbols0.value="ACDEFGHIKLMNPQRSTVWY";  document.form.color0.value="black";

            }
            else if (scheme=="amino_weblogo")  {
                document.form.symbols0.value="GSTYC";  document.form.color0.value="green";
                document.form.symbols1.value="NQ";  document.form.color1.value="purple";
                document.form.symbols2.value="KRH";  document.form.color2.value="blue";
                document.form.symbols3.value="DE";  document.form.color3.value="red";
                document.form.symbols4.value="PAWFLIMV";  document.form.color4.value="black";

            }
            else if (scheme=="amino_colors")  {
                document.form.symbols0.value="DE";  document.form.color0.value="bright red";
                document.form.symbols1.value="CM";  document.form.color1.value="yellow";
                document.form.symbols2.value="KR";  document.form.color2.value="blue";
                document.form.symbols3.value="ST";  document.form.color3.value="orange";
                document.form.symbols4.value="FY";  document.form.color4.value="mid blue";
                document.form.symbols5.value="NQ";  document.form.color5.value="cyan";
                document.form.symbols6.value="G";  document.form.color6.value="light grey";
                document.form.symbols7.value="LVI";  document.form.color7.value="green";
                document.form.symbols8.value="A";  document.form.color8.value="dark grey";
                document.form.symbols9.value="W";  document.form.color9.value="pink";
                document.form.symbols10.value="H"; document.form.color10.value="pale blue";
                document.form.symbols11.value="P"; document.form.color11.value="flesh";

            }
            else if (scheme=="amino_shapley")  {
                document.form.symbols0.value="DT";  document.form.color0.value="dark red";
                document.form.symbols1.value="E";  document.form.color1.value="red-brown";
                document.form.symbols2.value="C";  document.form.color2.value="bright yellow";
                document.form.symbols3.value="MY";  document.form.color3.value="dark yellow";
                document.form.symbols4.value="K";  document.form.color4.value="blue";
                document.form.symbols5.value="R";  document.form.color5.value="dark blue";
                document.form.symbols6.value="SQ";  document.form.color6.value="orange";
                document.form.symbols7.value="FPW";  document.form.color7.value="dark grey";
                document.form.symbols8.value="N";  document.form.color8.value="flesh";
                document.form.symbols9.value="GV";  document.form.color9.value="white";
                document.form.symbols10.value="I"; document.form.color10.value="dark green";
                document.form.symbols11.value="L"; document.form.color11.value="grey-green";
                document.form.symbols12.value="A"; document.form.color12.value="light green";
                document.form.symbols13.value="H"; document.form.color13.value="pale blue";

            }
            else if (scheme=="amino_charge")  {
                document.form.symbols0.value="KRH";  document.form.color0.value="blue";
                document.form.symbols1.value="DE";  document.form.color1.value="red";
                document.form.symbols2.value="ACFGILMNPQSTVWY";  document.form.color2.value="black";

            }
            else if (scheme=="amino_hydro")  {
                document.form.symbols0.value="AFGILPVWY";  document.form.color0.value="cyan";
                document.form.symbols1.value="CDEHKMNQRST";  document.form.color1.value="black";

            }
            else if (scheme=="amino_surface")  {
                document.form.symbols0.value="DEHKNPQRSTY";  document.form.color0.value="orange";
                document.form.symbols1.value="ACFGILMVW";  document.form.color1.value="black";

            }
            else if (scheme=="amino_flex")  {
                document.form.symbols0.value="DEKNPQRS";  document.form.color0.value="red";
                document.form.symbols1.value="ACFGHILMTVWY";  document.form.color1.value="green";

            }
            else if (scheme=="amino_disorder")  {
                document.form.symbols0.value="ARSQEGKP";  document.form.color0.value="red";
                document.form.symbols1.value="NCILFWYV";  document.form.color1.value="blue";
                document.form.symbols2.value="DHMT";  document.form.color2.value="black";

            }
            else if (scheme=="nucleo_bw")  {
                document.form.symbols0.value="ACGTU";  document.form.color0.value="black";

            }
            else if (scheme=="nucleo_weblogo")  {
                document.form.symbols0.value="G";  document.form.color0.value="orange";
                document.form.symbols1.value="TU";  document.form.color1.value="red";
                document.form.symbols2.value="C";  document.form.color2.value="blue";
                document.form.symbols3.value="A";  document.form.color3.value="green";

            }
            else if (scheme=="nucleo_shapley")  {
                document.form.symbols0.value="A";  document.form.color0.value="light blue";
                document.form.symbols1.value="C";  document.form.color1.value="orange";
                document.form.symbols2.value="G";  document.form.color2.value="light red";
                document.form.symbols3.value="T";  document.form.color3.value="light green";
                document.form.symbols4.value="U";  document.form.color4.value="dark grey";

            }
} 

