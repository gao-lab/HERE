# Template.rb - fills in the EPS template file with data, and
#               creates the EPS file. This file can be further
#               converted into PDF, GIF or PNG formats. 
#
# Vladimir Vacic, University of California, Riverside
# Lilia Iakoucheva, The Rockefeller University, New York
# Predrag Radivojac, Indiana University, Bloomington
#
# Oct-15-2005


class Template

    @@defaults = {
        'LOGO_HEIGHT'   => 5,
        'LOGO_WIDTH'    => 18,
        'LOGO_UNITS'    => 'cm',
        'YAXIS'         => 'false',
        'OUTLINE'       => 'false',
        'SHOW_BOX'      => 'n',
        'XAXIS'         => 'false',
        'SHRINK'        => 'false',
        'SHRINK_FACTOR' => 1,
        'FIXED_HEIGHT'  => 'false',
        'RES'           => 96,
        'RES_UNITS'     => 'ppi',
        'ANTIALIAS'     => true,
        'black'         => '[0 0 0]',
        'red'           => '[0.8 0 0]',
        'green'         => '[0 0.8 0]',
        'blue'          => '[0 0 0.8]',
        'yellow'        => '[1 0.71 0]',
        'purple'        => '[0.8 0 0.8]',
        'orange'        => '[1 0.7 0]'
    }


@@COLOR_SCHEME = {
    'amino_bw' => "
(G)  [0 0 0]
(S)  [0 0 0]
(T)  [0 0 0]
(Y)  [0 0 0]
(C)  [0 0 0]
(N)  [0 0 0]
(Q)  [0 0 0]
(K)  [0 0 0]
(R)  [0 0 0]
(H)  [0 0 0]
(D)  [0 0 0]
(E)  [0 0 0]
(P)  [0 0 0] 
(A)  [0 0 0]
(W)  [0 0 0]
(F)  [0 0 0]
(L)  [0 0 0]
(I)  [0 0 0]
(M)  [0 0 0]
(V)  [0 0 0]",

    'amino_weblogo' => "
(G)  [0 0.8 0] 
(S)  [0 0.8 0]
(T)  [0 0.8 0]
(Y)  [0 0.8 0]
(C)  [0 0.8 0]
(N)  [0.8 0 0.8] 
(Q)  [0.8 0 0.8]
(K)  [0 0 0.8] 
(R)  [0 0 0.8]
(H)  [0 0 0.8]
(D)  [0.8 0 0]
(E)  [0.8 0 0]
(P)  [0 0 0]
(A)  [0 0 0]
(W)  [0 0 0]
(F)  [0 0 0]
(L)  [0 0 0]
(I)  [0 0 0]
(M)  [0 0 0]
(V)  [0 0 0]",

    'amino_colors' => "
(G)  [0.9 0.9 0.9] 
(S)  [0.98 0.59 0]
(T)  [0.98 0.59 0]
(Y)  [0.2 0.2 0.67] 
(C)  [0.9 0.9 0] 
(N)  [0 0.86 0.86] 
(Q)  [0 0.86 0.86]
(K)  [0.08 0.35 1]
(R)  [0.08 0.35 1]
(H)  [0.51 0.51 0.82] 
(D)  [0.9 0.04 0.04] 
(E)  [0.9 0.04 0.04]
(P)  [0.86 0.59 0.71] 
(A)  [0.78 0.78 0.78] 
(W)  [0.71 0.35 0.71]  
(F)  [0.2 0.2 0.67]
(L)  [0.06 0.51 0.06] 
(I)  [0.06 0.51 0.06]
(M)  [0.9 0.9 0]
(V)  [0.06 0.51 0.06]",

    'amino_shapley' => "
(G)  [0.9 0.9 0.9]
(S)  [1 0.3 0.3]
(T)  [0.63 0 0.25]
(Y)  [0.72 0.63 0.26]
(C)  [1 1 0.44]
(N)  [1 0.49 0.44]
(Q)  [1 0.3 0.3]
(K)  [0.28 0.28 0.72]
(R)  [0 0 0.49]
(H)  [0.44 0.44 1]
(D)  [0.63 0 0.26]
(E)  [0.4 0 0]
(P)  [0.33 0.3 0.26]
(A)  [0.55 1 0.55]
(W)  [0.33 0.3 0.26]
(F)  [0.33 0.3 0.26]
(L)  [0.27 0.37 0.27]
(I)  [0 0.3 0]
(M)  [0.72 0.63 0.26]
(V)  [0.9 0.9 0.9]",

    'amino_charge' => "
(G)  [0 0 0]
(S)  [0 0 0]
(T)  [0 0 0]
(Y)  [0 0 0]
(C)  [0 0 0]
(N)  [0 0 0]
(Q)  [0 0 0]
(K)  [0 0 0.8] 
(R)  [0 0 0.8]
(H)  [0 0 0.8]
(D)  [0.8 0 0] 
(E)  [0.8 0 0]
(P)  [0 0 0]
(A)  [0 0 0]
(W)  [0 0 0]
(F)  [0 0 0]
(L)  [0 0 0]
(I)  [0 0 0]
(M)  [0 0 0]
(V)  [0 0 0]",

    'amino_hydro' => "
(G)  [0 1 1]
(S)  [0 0 0]
(T)  [0 0 0]
(Y)  [0 1 1]
(C)  [0 0 0]
(N)  [0 0 0]
(Q)  [0 0 0]
(K)  [0 0 0]
(R)  [0 0 0]
(H)  [0 0 0]
(D)  [0 0 0]
(E)  [0 0 0]
(P)  [0 1 1]
(A)  [0 1 1]
(W)  [0 1 1]
(F)  [0 1 1]
(L)  [0 1 1]
(I)  [0 1 1]
(M)  [0 0 0]
(V)  [0 1 1]",

    'amino_surface' => "
(G)  [0 0 0]
(S)  [1 0.65 0]
(T)  [1 0.65 0]
(Y)  [1 0.65 0]
(C)  [0 0 0]
(N)  [1 0.65 0]
(Q)  [1 0.65 0]
(K)  [1 0.65 0]
(R)  [1 0.65 0]
(H)  [1 0.65 0]
(D)  [1 0.65 0]
(E)  [1 0.65 0]
(P)  [1 0.65 0]
(A)  [0 0 0]
(W)  [0 0 0]
(F)  [0 0 0]
(L)  [0 0 0]
(I)  [0 0 0]
(M)  [0 0 0]
(V)  [0 0 0]",

    'amino_flex' => "
(G)  [0 0.8 0]
(S)  [0.8 0 0]
(T)  [0 0.8 0]
(Y)  [0 0.8 0]
(C)  [0 0.8 0]
(N)  [0.8 0 0]
(Q)  [0.8 0 0]
(K)  [0.8 0 0]
(R)  [0.8 0 0]
(H)  [0 0.8 0]
(D)  [0.8 0 0]
(E)  [0.8 0 0]
(P)  [0.8 0 0]
(A)  [0 0.8 0]
(W)  [0 0.8 0]
(F)  [0 0.8 0]
(L)  [0 0.8 0]
(I)  [0 0.8 0]
(M)  [0 0.8 0]
(V)  [0 0.8 0]",

    'amino_disorder' => "
(G)  [0.8 0 0]
(S)  [0.8 0 0]
(T)  [0 0 0]
(Y)  [0 0 0.8]
(C)  [0 0 0.8]
(N)  [0 0 0.8]
(Q)  [0.8 0 0]
(K)  [0.8 0 0]
(R)  [0.8 0 0]
(H)  [0 0 0]
(D)  [0 0 0]
(E)  [0.8 0 0]
(P)  [0.8 0 0]
(A)  [0.8 0 0]
(W)  [0 0 0.8]
(F)  [0 0 0.8]
(L)  [0 0 0.8]
(I)  [0 0 0.8]
(M)  [0 0 0]
(V)  [0 0 0.8]",

    'nucleo_bw' => "
(G)  [0 0 0]
(T)  [0 0 0]
(C)  [0 0 0]
(A)  [0 0 0]
(U)  [0 0 0]",

    'nucleo_weblogo' => "
(G)  [1 0.7 0] 
(T)  [0.8 0 0]
(C)  [0 0 0.8] 
(A)  [0 0.8 0]
(U)  [0.8 0 0]",

    'nucleo_shapley' => "
(G)  [1 0.44 0.44]
(T)  [0.63 1 0.63]
(C)  [1 0.55 0.29]
(A)  [0.63 0.63 1]
(U)  [0.72 0.72 0.72]"
}

    # Creates the EPS file accoring to the data and a set of parameters,
    # and if the output format is different from EPS, converts it into
    # GIF, PNG, or PDF.
    # 
    # Input: conserved - array of conserved symbols 
    #        enriched  - array of heights for enriched symbols
    #        depleted  - array of heights for depleted symbols
    #        options   - hash of parameters passed from the driver script 
    #        outfile   - the name of the output file
    #        path      - output directory (write permissions should be set 
    #                    accordingly) 
    def Template.fill_template(conserved, enriched, depleted, options, outfile, path)

        # put parameters in params hash table
        params = make_params(options)

        # set default data if not filled
        params = set_defaults(params, enriched.length)

        # put color in params
        params = set_colors(params)

        # put data in params
        params = set_data(params, conserved, enriched, depleted)

        # make EPS output
        eps    = fill("template.eps", params)
        format = options['FORMAT']

        # convert
        answer  = get_configuration
        gs      = answer[0]
        convert = answer[1]
     
        width     = params['LOGO_WIDTH_POINTS']
        height    = params['LOGO_HEIGHT_POINTS']    
        res       = params['RES']

        if params['RES_UNITS']=="ppc"
            res *= 2.54
        elsif params['RES_UNITS']=="ppp"
            res *= 72 
        end

        antialias = (options['ANTIALIAS'] != nil && options['ANTIALIAS']) ? "-dTextAlphaBits=4" : ""

        r = "#{width}x#{height}"

        if format == "EPS"
            if outfile == "-"  # if standard out
                puts eps
            else
                file = File.open("#{path}#{outfile}", "w")
                file.puts(eps)
                file.close
            end

        elsif format == "PDF"
            program = "#{gs} -sOutputFile=#{path}#{outfile} -sDEVICE=pdfwrite -dPDFSETTINGS=/printer -r#{res} -dDEVICEWIDTHPOINTS=#{width} -dDEVICEHEIGHTPOINTS=#{height} -dEmbedAllFonts=true -dSAFER -dBATCH -dNOPAUSE -_"
            f = IO.popen(program, "w")
            f.puts(eps)
            f.close

        elsif format == "PNG"
            program = "#{gs} -sOutputFile=#{path}#{outfile} -sDEVICE=png16m -q -r#{res} -dDEVICEWIDTHPOINTS=#{width} -dDEVICEHEIGHTPOINTS=#{height} #{antialias} -dSAFER -dBATCH -dNOPAUSE -_"
            f = IO.popen(program, "w+")
            f.puts(eps)
            f.close

        elsif format == "GIF"  # convert to PNG first, then GIF
            program = "#{gs} -sOutputFile=- -sDEVICE=png16m -q -r#{res} -dDEVICEWIDTHPOINTS=#{width} -dDEVICEHEIGHTPOINTS=#{height} #{antialias} -dSAFER -dBATCH -dNOPAUSE -_"

            if outfile == "-"
                program += " | #{convert} png:- gif:-"
            else
                program += " | #{convert} png:- #{path}#{outfile}"
            end

            f = IO.popen(program, "w+")
            f.puts(eps)
            f.close
        end
    end


    #
    # Reads the configuration file (tsl.conf).
    #
    # Output: answer[0] - path to GhostScript binary gs
    #         answer[1] - path to ImageMagick binary convert
    #
    def Template.get_configuration
        answer = Array.new
        answer << "gs" << "convert"     

        # No configuration file, then return defaults.
        return answer if !FileTest.exists?("tsl.conf")

        config = File.open("tsl.conf")
    
        config.each do |line|
            next if line =~ /^\#/       # skip lines beginning with "#"

            if line =~ /^gs/i           # if looks like gs (case insensitive)
                line =~ /^\S+\=(.+)$/
                answer[0] = $1
            elsif line =~ /^convert/i   # if looks like convert (case insensitive)
                line =~ /^\S+\=(.+)$/
                answer[1] = $1
            end
        end

        if answer[0] == nil || !FileTest.exists?(answer[0])
            puts "Please check tsl.conf: gs program (#{answer[0]}) does not exist."
            exit
        end

        if answer[1] == nil || !FileTest.exists?(answer[1])
            puts "Please check tsl.conf: convert program (#{answer[1]}) does not exist." 
            exit
        end

        return answer
    end


    #
    # Fills in the EPS template with acutal values.
    #
    # Input:  filename - EPS template file name
    #         params   - hash table with template substitutions
    #
    # Output: text     - array of lines after substitution
    #
    def Template.fill(filename, params)
        if !FileTest.exists?(filename)
	    puts "Template filename (#{filename}) does not exist.\n"
            exit
        end

        file = File.open(filename)    
        text = file.readlines
        file.close

        #replace {$KEYWORDS} with value in %$params hash
        text.each_index do |i|
            text[i] = text[i].gsub(/\{\$(.*?)\}/)  {  |match|
                params[$1] != nil ? params[$1] : ""
            }
        end

        return text 
    end


    #
    # Puts user parameters into a local hash.
    #
    def Template.make_params(options)
        params = {
            'INPUT_KIND'    => options['INPUT_KIND'],
            'LOGO_HEIGHT'   => options['LOGO_HEIGHT'].to_f,
            'LOGO_WIDTH'    => options['LOGO_WIDTH'].to_f,
            'LOGO_UNITS'    => options['LOGO_UNITS'],
            'OUTLINE'       => options['OUTLINE'],
            'XAXIS'         => options['XAXIS'],
            'LOGO_START'    => options['LOGO_START'],
            'LOGO_END'      => options['LOGO_END'],
            'START_NUM'     => options['START_NUM'],
            'YAXIS'         => options['YAXIS'],
            'TITLE'         => options['TITLE'],
            'SHOW_BOX'      => options['SHOW_BOX'] ? "s" : "n",
            'SHRINK'        => options['SHOW_BOX'],
            'SHRINK_FACTOR' => options['SHRINK_FACTOR'],
            'MAX_HEIGHT'    => options['MAX_HEIGHT'],
            'FIXED_HEIGHT'  => options['FIXED_HEIGHT'],
            'COLOR_SCHEME'  => options['COLOR_SCHEME'],
            'USER_COLORS'   => options['USER_COLORS'],
            'RES'           => options['RES'],
            'RES_UNITS'     => options['RES_UNITS']
        }
    
        params['MAX_LABEL'] = ("%4.1f" % (options['MAX_HEIGHT']*100)) << "%"

        return params
    end


    #
    # Fills in the missing parameters with the default values.
    #
    def Template.set_defaults(params, numchars)
        params['CHARS_PER_LINE'] = numchars

        params['LOGO_HEIGHT'] = @@defaults['LOGO_HEIGHT'] if params['LOGO_HEIGHT']==0 
        params['LOGO_WIDTH']  = @@defaults['LOGO_WIDTH']  if params['LOGO_WIDTH']==0

        params['RES']       = @@defaults['RES']       if params['RES'] == nil
        params['RES_UNITS'] = @@defaults['RES_UNITS'] if params['RES_UNITS'] == nil
        params['ANTIALIAS'] = @@defaults['ANTIALIAS'] if params['ANTIALIAS'] == nil
        
        if params['LOGO_UNITS']=="inch"
            params['LOGO_HEIGHT'] *= 2.54
            params['LOGO_WIDTH']  *= 2.54
        elsif params['LOGO_UNITS']=="pixel"
            params['LOGO_HEIGHT'] *= 2.54 / params['RES']
            params['LOGO_WIDTH']  *= 2.54 / params['RES']
        elsif params['LOGO_UNITS']=="point"
            params['LOGO_HEIGHT'] *= 2.54 / 72
            params['LOGO_WIDTH']  *= 2.54 / 72
        end

        if params['LOGO_START']==nil || params['LOGO_START']==""
            params['LOGO_START'] = params['START_NUM']
        else
            params['LOGO_START'] = params['LOGO_START'].to_i
        end

        if params['LOGO_END']==nil || params['LOGO_END'] == ""
            params['LOGO_END'] = numchars + params['LOGO_START'] - 1
        else
            params['LOGO_END'] = params['LOGO_END'].to_i
        end
  
        params['TITLE'] = "" if params['TITLE'] == nil

        params['YAXIS'] = @@defaults['YAXIS'] if params['YAXIS'] == nil
        params['XAXIS'] = @@defaults['XAXIS'] if params['XAXIS'] == nil

        params['OUTLINE']      = @@defaults['OUTLINE'] if params['OUTLINE'] == nil
        params['SHOW_BOX']     = @@defaults['SHOW_BOX'] if params['SHOW_BOX'] == nil
        params['FIXED_HEIGHT'] = @@defaults['FIXED_HEIGHT'] if params['FIXED_HEIGHT'] == nil

        params['SHRINK'] = @@defaults['SHRINK'] if params['SHRINK_FACTOR'] == nil
        params['SHRINK_FACTOR'] = @@defaults['SHRINK_FACTOR'] if params['SHRINK_FACTOR'] == nil

        givenrange    = params['LOGO_END'] - params['LOGO_START'] + 1
        possiblerange = numchars - (params['LOGO_START'] - params['START_NUM'])

        params['CHAR_WIDTH'] = (params['LOGO_WIDTH'] - 1.5) / numchars 

        params['LOGO_HEIGHT_POINTS'] = (params['LOGO_HEIGHT'] * (72 / 2.54)).round
        params['LOGO_WIDTH_POINTS']  = (params['LOGO_WIDTH']  * (72 / 2.54)).round
        params['LOGO_LINE_HEIGHT']   = params['LOGO_HEIGHT'] / 2.0

        return params
    end


    #
    # Adds residue color definitions to the hash of paramters.
    #
    def Template.set_colors(params)
        color_dict = "/colorDict <<\n"

        if params['COLOR_SCHEME'] != "user_colors"
            params['COLOR_SCHEME'] = "amino_weblogo" if params['COLOR_SCHEME']==nil || params['COLOR_SCHEME']==""
            color_dict += @@COLOR_SCHEME[params['COLOR_SCHEME']]

        else
            params['USER_COLORS'].keys.each do |key|
                color_dict += "(#{key}) #{params['USER_COLORS'][key]}\n"
            end

        end

        color_dict += "\n>> def"

        params['COLOR_DICT'] = color_dict

        return params
    end


    #
    # Sets the DATA field in the EPS template.
    # 
    # Input:  params         - hash table of parameters
    #         conserved_data -
    #         enriched_data  -
    #         depleted_data  -
    #
    # Output: params - hash table of parameters 
    #
    def Template.set_data(params, conserved_data, enriched_data, depleted_data)
        height = ""
        letter = ""

        conserved = ""
        enriched  = "" 
        depleted  = ""

        start_num  = params['START_NUM']
        logo_start = params['LOGO_START'] - start_num    # where in data to start
        logo_end   = params['LOGO_END']   - start_num    # where in data to end
        numlabel   = params['LOGO_START']

        logo_end = (logo_end >= enriched_data.length) ? (enriched_data.length - 1) : logo_end;

        (logo_start..logo_end).each do |i|
            #
            # Conserved residues
            #
            conserved += <<END
(#{numlabel}) StartStack
END

            conserved += " middleMargin stackMargin sub pointsPerBit div (#{conserved_data[i]}) MakeSymbol\n"

            conserved += <<END
EndStack

END

            #
            # Enriched residues
            #
            slice = enriched_data[i]
            enriched += <<END
(#{numlabel}) StartStack
END
    

            slice.each do |s|
                s =~ /^(.{1})(\S+)/
                letter = $1.upcase
                height = $2.to_f

                # is space, so leave
                break if letter == " "
	    
                enriched += " #{height} (#{letter}) MakeSymbol\n"
            end

            enriched += <<END
EndStack

END

            #
            # Depleted residues
            #
            slice = depleted_data[i]
            depleted += <<END
(#{numlabel}) StartStack
END

            slice.each do |s|
                s =~ /^(.{1})(\S+)/
                letter = $1.upcase
                height = $2.to_f

                # is space, so leave
                break if letter == " "

                depleted += " #{height} (#{letter}) MakeSymbol\n"
            end

            depleted += <<END
EndStack

END

            numlabel += 1
        end

        params['DATA_CONSERVED'] = conserved
        params['DATA_ENRICHED']  = enriched
        params['DATA_DEPLETED']  = depleted
   
        return params
    end

end
