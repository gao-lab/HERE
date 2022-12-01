#!/usr/bin/env ruby

# CGI interface to the TSL classes.
#
# Version 1.2
#
# Vladimir Vacic, University of California, Riverside
# Lilia M. Iakoucheva, The Rockefeller University, New York
# Predrag Radivojac, Indiana University, Bloomington
# 
# Nov-27-2006

require "cgi"
require "TSL.rb"
require "Template.rb"


path = "../../"
temp = "/var/www/html/tsl/cache/"

cgi = CGI.new("html3")

command = ""
command = cgi["command"].read if cgi.has_key?("command")

positive_sample = ""
positive_sample = cgi["positive_sample"].read if cgi.has_key?("positive_sample")

positive_file = ""
positive_file = cgi["positive_file"].read if cgi.has_key?("positive_file")

negative_sample = ""
negative_sample = cgi["negative_sample"].read if cgi.has_key?("negative_sample")

negative_file = ""
negative_file = cgi["negative_file"].read if cgi.has_key?("negative_file")

input_kind = "A"
input_kind = cgi["input_kind"].read if cgi.has_key?("input_kind")

test = "ttest"
test = cgi["test"].read if cgi.has_key?("test")

p_value = 0.05
p_value = cgi["p_value"].read.to_f if cgi.has_key?("p_value")

bonferroni = false 
bonferroni = cgi["bonferroni"].read=="on" if cgi.has_key?("bonferroni") 

title = ""
title = cgi["title"].read if cgi.has_key?("title")

format = "PNG"
format = cgi["format"].read if cgi.has_key?("format")

fixed_height = false 
fixed_height = ("on" == cgi["fixed_height"].read) if cgi.has_key?("fixed_height") 

logo_height = "5" 
logo_height = cgi["logo_height"].read if cgi.has_key?("logo_height")

logo_width = "18"
logo_width = cgi["logo_width"].read if cgi.has_key?("logo_width")

logo_units = "cm"
logo_units = cgi["logo_units"].read if cgi.has_key?("logo_units")

start_num = 1 
start_num = cgi["start_num"].read.to_i if cgi.has_key?("start_num")

logo_start = "" 
logo_start = cgi["logo_start"].read if cgi.has_key?("logo_start")

logo_end = ""
logo_end = cgi["logo_end"].read if cgi.has_key?("logo_end")

xaxis = false
xaxis = cgi["xaxis"].read=="on" if cgi.has_key?("xaxis")

yaxis = false
yaxis = cgi["yaxis"].read=="on" if cgi.has_key?("yaxis")

box = false
box = cgi["box"].read=="on" if cgi.has_key?("box")

shrink = ""
shrink= cgi["shrink"].read if cgi.has_key?("shrink")

antialias = false
antialias = ("on"==cgi["antialias"].read)  if cgi.has_key?("antialias")

res = 96
res = cgi["res"].read.to_i if cgi.has_key?("res")

res_units = "ppi"
res_units = cgi["res_units"].read if cgi.has_key?("res_units")

outline = false
outline = cgi["outline"].read=="on" if cgi.has_key?("outline")

colors = "amino_weblogo"
colors = "nucleo_weblogo"  if ("N" == input_kind)
colors = cgi["colors"].read if cgi.has_key?("colors")


$header=<<END
<html>
    <head>
        <title>Two Sample Logos - Create</title>
        <link rel="stylesheet" href="#{path}tsl.css" type="text/css">
        <link rel="shortcut icon" href="#{path}images/favicon.ico">
        <meta http-equiv="Content-Type" content="text/html;charset=UTF-8">
        <script language="JavaScript" type="text/javascript" src="#{path}set_colors.js"></script>
    </head>
    <body>
        <table width="100%" cellpadding="0" cellspacing="6" border="0" align="center">
        <tr>
        <td width="50%" align="left" valing="bottom"><a href="#{path}index.html"><img src="#{path}images/tsl.gif" width="206" height="60" border="0" alt="Two Sample Logos" hspace="0" vspace="6"></a></td>
        <td align="right" valign="bottom">
<a href="#{path}index.html">Home</a> |
Create TSL| 
<a href="#{path}examples.html">Examples</a> |
<a href="#{path}help.html">Help</a></td>
        </tr>
END

$body=<<END
        <form name="form" method="post" action="tsl.cgi" enctype="multipart/form-data">

        <tr>
        <td width="50%" valign="top" align="center">
        <h5>Positive Sample</h5>

            <table><tr><td>
            <textarea name="positive_sample" rows="8" cols="55">#{positive_sample}</textarea><br>
            Upload positive sample file: <input type="file" name="positive_file">
            <a href="#{path}help.html#sequence"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a><br><br>
            </td></tr></table>

        </td>
        <td width="50%" valign="top" align="center">
        <h6>Negative Sample</h6>

            <table><tr><td>
            <textarea name="negative_sample" rows="8" cols="55">#{negative_sample}</textarea><br>
            Upload negative sample file: <input type="file" name="negative_file">
            <a href="#{path}help.html#sequence"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a><br><br>
            </td></tr></table>

        </td>
        </tr>

        <tr>
        <td colspan="2">
            <table width="100%" cellspacing="2" cellpadding="0">
            <tr>
            <td width="33%" valign="top">
            <h4>Basic Options</h4>
            </td>
            <td rowspan="8">&nbsp; </td>
            <td width="33%" valign="top">
            <h4>Advanced Options</h4>
            </td>
            <td rowspan="8">&nbsp; </td>
            <td width="33%" valign="top">
            <h4>Output Options</h4>
            </td>
            </tr>


            <tr>
            <td width="33%">&nbsp; Sequence Type:
            <input type="radio" name="input_kind" value="A" #{"checked=\"checked\"" if input_kind=="A"}> AA 
            <input type="radio" name="input_kind" value="N" #{"checked=\"checked\"" if input_kind=="N"}> DNA / RNA 
            <a href="#{path}help.html#seq_type"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a>
            </td>
            <td width="33%">&nbsp; Title:
            <input type="text" name="title" size="20" maxlength="80" value="#{title}">
            <a href="#{path}help.html#title"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a>
            </td>
            <td width="33%">&nbsp; Output format: 
            <select name="format">
                <option value="GIF">GIF (bitmap)</option>
                <option value="PNG" selected="selected">PNG (bitmap)</option>
                <option value="EPS">EPS (vector)</option>
                <option value="PDF">PDF (vector)</option>
                <option value="TXT">TXT (raw values)</option>
            </select>
            <a href="#{path}help.html#format"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a>
            </td>
            </tr>


            <tr>
            <td width="33%">&nbsp; Statistical test:
            <select name="test">
                <option value="ttest" selected="selected">t-test</option>
                <option value="binomial">binomial test</option>
            </select>
            <a href="#{path}help.html#test"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a>
            </td>
            <td width="33%">&nbsp; Logo range:
            <input type="text" name="logo_start" size="4" maxlength="80" value="#{logo_start}"> -
            <input type="text" name="logo_end" size="4" maxlength="80" value="#{logo_end}">
            <a href="#{path}help.html#range"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a>
            </td>
            <td width="33%">&nbsp; Output size:
            <input type="text" name="logo_width" value="18" size="6" maxlength="80" value="#{logo_width}"> X
            <input type="text" name="logo_height" value="5" size="6" maxlength="80" value="#{logo_height}"> 
            <select name="logo_units">
                <option value="cm" selected="selected">cm</option>
                <option value="inches">inches</option>
                <option value="pixels">pixels</option>
                <option value="points">points</option>
            </select>
            <a href="#{path}help.html#size"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a>
            </td>
            </tr>


            <tr>
            <td width="33%">&nbsp; Show residue if p-value &lt;
            <select name="p_value">
                <option value="0.0001">0.0001</option>
                <option value="0.005">0.005</option>
                <option value="0.01">0.01</option>
                <option value="0.05" selected="selected">0.05</option>
                <option value="0.1">0.1</option>
                <option value="0.5">0.5</option>
                <option value="1">1</option>
            </select>
            <a href="#{path}help.html#p_value"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a>
            </td>
            <td>&nbsp; First position index:
            <input type="text" name="start_num" value="#{start_num}" size="4" maxlength="80">
            <a href="#{path}help.html#firstnum"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a>
            </td>
            <td width="33%">&nbsp; Resolution:
            <input type="text" name="res" value="96" size="4" maxlength="8">
            <select name="res_units">
                <option value="ppc">pixels/cm</option>
                <option selected="selected" value="ppi">pixels/inch (dpi)</option>
                <option value="ppp">pixels/point</option>
            </select>
            <a href="#{path}help.html#resolution"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a>
            </td>
            </tr>


            <tr>
            <td>&nbsp; Bonferroni correction:
            <input type="checkbox" name="bonferroni" value="on">
            <a href="#{path}help.html#bonferroni"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a>
            </td>
            <td width="33%">&nbsp; Show X-axis indexes:
            <input type="checkbox" name="xaxis" value="on" checked="checked">
            <a href="#{path}help.html#xaxis"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a>
            </td>
            <td>&nbsp; Antialias bitmaps:
            <input type="checkbox" name="antialias" value="on" checked="checked">
            <a href="#{path}help.html#antialias"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a>
            </td>
            </tr>


            <tr>
            <td width="33%">&nbsp; Show conserved residues:
            <input type="checkbox" name="conserved" value="on" checked="checked">
            <a href="#{path}help.html#conserved"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a>
            </td>
            <td>&nbsp; Show Y-axis labels:
            <input type="checkbox" name="yaxis" value="on" checked="checked">
            <a href="#{path}help.html#yaxis"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a>
            </td>
            <td>&nbsp; Boxed / Boxed shrink factor:
            <input type="checkbox" name="box" value="on"> 
            <input type="text" name="shrink" value="0.5" size="4" maxlength="80" /> 
            <a href="#{path}help.html#box"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a>
            </td>
            </tr> 


            <tr>
            <td width="33%">&nbsp; Fixed height symbols:
            <input type="checkbox" name="fixed_height" value="on">
            <a href="#{path}help.html#fixed"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a>
            </td>
            <td>&nbsp; 
            </td>
            <td>&nbsp; Outline symbols:
            <input type="checkbox" name="outline" value="on">
            <a href="#{path}help.html#outline"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a>
            </td>
            </tr> 


            <tr>
            <td width="33%" valign="top">
            <hr noshade size="1">
            </td>
            <td width="33%" valign="top">
            <hr noshade size="1">
            <center>
            <input type="hidden" name="command" value="create">
            <input type="submit" name="submit"  value="Create TSL">
            <input type="reset"  name="reset"   value="Reset">
            </center>
            </td>
            <td width="33%" valign="top">
            <hr noshade size="1">
            </td>
            </tr>
            </table>
        </td>
        </tr>

        <tr>
        <td colspan="2">
            <h4>Color Scheme</h4>
        </td>
        </tr>

        <tr>
        <td colspan="2">

            <table width="100%">

            <tr>
            <td width="33%" valign="top"><b>Amino acids:</b><br>
            <input type="radio" name="colors" value="amino_bw" onClick="setColors(this)">Black and white <a href="#{path}help.html#bw"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a><br>
            <input type="radio" name="colors" value="amino_weblogo" #{"checked=\"checked\"" if input_kind=="A"} onClick="setColors(this)">WebLogo defaults <a href="#{path}help.html#amino_weblogo"><img src="#{path}images/info.gif" width="12" height="12" border="0"></a><br>
            <input type="radio" name="colors" value="amino_colors" onClick="setColors(this)">Amino Color <a href="#{path}help.html#amino_colors"><img src="#{path}images/info.gif" width="12" height="12" border="0"></a><br>
            <input type="radio" name="colors" value="amino_shapley" onClick="setColors(this)">Shapley <a href="#{path}help.html#amino_shapley"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a><br>
            </td>

            <td width="33%" valign="top"><b>Amino acid properties:</b><br>
            <input type="radio" name="colors" value="amino_charge" onClick="setColors(this)">Charge <a href="#{path}help.html#amino_charge"><img src="#{path}images/info.gif" width="12" height="12" border="0"></a><br>
            <input type="radio" name="colors" value="amino_hydro" onClick="setColors(this)">Hydrophobicity <a href="#{path}help.html#amino_hydro"><img src="#{path}images/info.gif" width="12" height="12" border="0"></a><br>
            <input type="radio" name="colors" value="amino_surface" onClick="setColors(this)">Surface exposure <a href="#{path}help.html#amino_surface"><img src="#{path}images/info.gif" width="12" height="12" border="0"></a><br>
            <input type="radio" name="colors" value="amino_flex" onClick="setColors(this)">Flexibility <a href="#{path}help.html#amino_flex"><img src="#{path}images/info.gif" width="12" height="12" border="0"></a><br>
            <input type="radio" name="colors" value="amino_disorder" onClick="setColors(this)">Disorder <a href="#{path}help.html#amino_disorder"><img src="#{path}images/info.gif" width="12" height="12" border="0"></a><br>
            </td>

            <td width="33%" valign="top"><b>Nucleotides:</b><br>
            <input type="radio" name="colors" value="nucleo_bw" onClick="setColors(this)">Black and white <a href="#{path}help.html#bw"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a><br>
            <input type="radio" name="colors" value="nucleo_weblogo" #{"checked=\"checked\"" if input_kind=="N"} onClick="setColors(this)">WebLogo defaults <a href="#{path}help.html#nucleo_weblogo"><img src="#{path}images/info.gif" width="12" height="12" border="0"></a><br>
            <input type="radio" name="colors" value="nucleo_shapley" onClick="setColors(this)">Shapley <a href="#{path}help.html#nucleo_shapley"><img src="#{path}images/info.gif" width="12" height="12" border="0" alt="help"></a><br>

            <br><b>Other:</b><br>
            <input type="radio" name="colors" value="user_colors" onClick="setColors(this)">User defined <a href="#{path}help.html#user_colors"><img src="#{path}images/info.gif" width="12" height="12" border="0"></a><br>
            </td>

            </tr>
            </table>

            <br><hr noshade size="1"><br>

            <table width="100%">
            <tr>
            <td width="16%"><b>Symbols:</b></td>
            <td width="17%"><b>Color:</b></td>
            <td width="16%"><b>Symbols:</b></td>
            <td width="17%"><b>Color:</b></td>
            <td width="16%"><b>Symbols:</b></td>
            <td width="17%"><b>Color:</b></td>
            </tr>

            <tr>
            <td width="16%" valign="top">
            <input type="text" name="symbols0"><br>
            <input type="text" name="symbols1"><br>
            <input type="text" name="symbols2"><br>
            <input type="text" name="symbols3"><br>
            <input type="text" name="symbols4"><br>
            </td>

            <td width="17%" valign="top">
            <input type="text" name="color0"><br>
            <input type="text" name="color1"><br>
            <input type="text" name="color2"><br>
            <input type="text" name="color3"><br>
            <input type="text" name="color4"><br>
            </td>

            <td width="16%" valign="top">
            <input type="text" name="symbols5"><br>
            <input type="text" name="symbols6"><br>
            <input type="text" name="symbols7"><br>
            <input type="text" name="symbols8"><br>
            <input type="text" name="symbols9"><br>
            </td>

            <td width="17%" valign="top">

            <input type="text" name="color5"><br>
            <input type="text" name="color6"><br>
            <input type="text" name="color7"><br>
            <input type="text" name="color8"><br>
            <input type="text" name="color9"><br>
            </td>

            <td width="16%" valign="top">
            <input type="text" name="symbols10"><br>
            <input type="text" name="symbols11"><br>
            <input type="text" name="symbols12"><br>
            <input type="text" name="symbols13"><br>
            <input type="text" name="symbols14"><br>
            </td>

            <td width="17%" valign="top">
            <input type="text" name="color10"><br>
            <input type="text" name="color11"><br>
            <input type="text" name="color12"><br>
            <input type="text" name="color13"><br>
            <input type="text" name="color14"><br>
            </td>
            </tr>
            </table>

        </td>
        </tr>

        <tr>
        <td colspan="2">
            <hr noshade size="1">
            <center>
            <input type="hidden" name="command" value="create">
            <input type="submit" name="submit"  value="Create TSL">
            <input type="reset"  name="reset"   value="Reset">
            </center>
            <hr noshade size="3">
        </td>
        </tr>

        </form>
END


$footer=<<END
        </table>
    </body>
</html>
END


#
# Utility function for displaying error messages.
#
def show_error(message)
    print "Content-type: text/html\n\n"
    print $header
    print "<tr><td colspan=\"2\" align=\"center\"><br><br><br><font color=\"red\" size=\"+1\">Error: #{message}</font></td></tr>\n"  
    print $footer
    exit
end


#
# Main code starts here.
#

case command
when "create"
    pos = Array.new
    if (positive_sample.empty?)
        if (positive_file.empty?)
            show_error("Positive sample missing.")
        else
            pos = positive_file.upcase.split("\n")
        end
    else
        pos = positive_sample.upcase.split("\n")
    end


    neg = Array.new
    if (negative_sample.empty?)
        if (negative_file.empty?)
            show_error("Negative sample missing.")
        else
            neg = negative_file.upcase.split("\n")
        end
    else
        neg = negative_sample.upcase.split("\n")
    end


    ########################################################################
    #
    # Read input files and compute enrichment/depletion and p-values 
    #
    ########################################################################

    params = {
        'INPUT_KIND'   => input_kind,
        'TEST'         => test,
        'P_VALUE'      => p_value,
        'BONFERRONI'   => bonferroni,
        'FIXED_HEIGHT' => fixed_height
    }

    result = TSL.get_p_values(pos, neg, params)

    error_code  = result[0]
    error_line  = result[1]
    seq_length  = result[2]
    p_values    = result[3]
    frequencies = result[4]


    ########################################################################
    #
    # Write output, if there were no errors 
    #
    ########################################################################

    case error_code
    when 1
        show_error("Positive example #{error_line} is malformatted.")

    when 2
        show_error("Positive example #{error_line} differs in length.")

    when 3
        show_error("Negative example #{error_line} is malformatted.")

    when 4
        show_error("Negative example #{error_line} differs in length.")

    when 5
        show_error("Positive and negative samples differ in sequence length.")

    when 0
        if ("TXT" == format)
            print "Content-type: text/plain\n\n"

            print "#POSITION\tSYMBOL\tFREQUENCY\tP-VALUE\n" 

            for i in 0...(seq_length)
                if ("A" == input_kind)
                    for j in 0...(TSL.AMINOACIDS.length)
                        if (p_values[i][j] < p_value)
                            print "#{i + start_num}\t#{TSL.AMINOACIDS[j]}\t#{frequencies[i][j]}\t#{p_values[i][j]}\n"
                        end
                    end
                else
                    for j in 0...(TSL.NUCLEOTIDES.length)
                        if (p_values[i][j] < p_value)
                            print "#{i + start_num}\t#{TSL.NUCLEOTIDES[j]}\t#{frequencies[i][j]}\t#{p_values[i][j]}\n"
                        end
                    end
                end
            end
        else
            result = TSL.get_heights(seq_length, p_values, frequencies, params)

            conserved   = result[0]
            enriched    = result[1]
            depleted    = result[2]
            max_height  = result[3]

            max_height = 0.05  if (0==max_height)

            user_colors = Hash.new

            if ("user_colors" == colors)
                (0..14).each do |i|
                    symbols = cgi["symbols#{i}"].read if cgi.has_key?("symbols#{i}")
                    color = cgi["color#{i}"].read if cgi.has_key?("color#{i}")

                    ("A".."Z").each do |letter|
                        user_colors[letter] = color if symbols.include?(letter)
                    end
                end
            end

            options = {
                'INPUT_KIND'    => input_kind,
                'LOGO_HEIGHT'   => logo_height,
                'LOGO_WIDTH'    => logo_width,
                'LOGO_UNITS'    => logo_units,
                'COLOR_SCHEME'  => colors,
                'USER_COLORS'   => user_colors, 
                'LOGO_START'    => logo_start,
                'LOGO_END'      => logo_end,
                'START_NUM'     => start_num,
                'TITLE'         => title,
                'XAXIS'         => xaxis,
                'YAXIS'         => yaxis,
                'SHRINK_FACTOR' => shrink,
                'FIXED_HEIGHT'  => fixed_height,
                'MAX_HEIGHT'    => max_height,
                'RES'           => res,
                'RES_UNITS'     => res_units,
                'FORMAT'        => format,
                'ANTIALIAS'     => antialias,
	        'OUTLINE'       => outline,
                'SHOW_BOX'      => box
            }

            output = "#{rand.to_s[2,15].to_s}.#{format.downcase}"

            Template.fill_template(conserved, enriched, depleted, options, output, temp)

            case format
            when "PNG"
                # print "Content-Type: image/png\n\n"
                print "Content-type: text/html\n\n"
                print "<html><body><img src=\"../../cache/#{output}\"></body></html>"
            when "GIF"
                # print "Content-Type: image/gif\n\n"
                print "Content-type: text/html\n\n"
                print "<html><body><img src=\"../../cache/#{output}\"></body></html>"
            when "EPS"
                print "Content-Type: application/eps\n\n"
                print `cat #{temp}#{output}`
            when "PDF"
                print "Content-Type: application/pdf\n\n"
                print `cat #{temp}#{output}`
            end
        end
    end
else
    print "Content-type: text/html\n\n"
    print $header
    print $body
    print $footer
end
