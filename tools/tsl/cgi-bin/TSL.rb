# Processes sequences from FASTA, CLUSTAL or flat files and computes enrichment/
# depletion and p-values.
#
# Version 1.2
#
# Vladimir Vacic, University of California, Riverside
# Lilia Iakoucheva, The Rockefeller University, New York
# Predrag Radivojac, Indiana University, Bloomington
#
# Oct-15-2005, May-1-2009


class TSL 
    @@NUCLEOTIDES = %w{A G C T U}
    @@AMINOACIDS  = %w{A C D E F G H I K L M N P Q R S T V W Y}

    def TSL.NUCLEOTIDES
        @@NUCLEOTIDES
    end

    def TSL.AMINOACIDS
        @@AMINOACIDS
    end


    # Input:  pos_data    - array of lines from the positive sample
    #         neg_data    - array of lines from the negative sample 
    #         params      - hash table of parameters
    #
    # Output: An array containing the following:
    #         answer[0]   - 0    everything is OK
    #                       1    positives file is malformatted
    #                       2    not all positive sequences have the same length
    #                       3    negatives file is malformatted
    #                       4    not all negative sequences have the same length
    #                       5    positives and negatives do not have the same length
    #         answer[1]   - location of the error, or nil if error code is 0
    #         answer[2]   - sequence length
    #         answer[3]   - a sequence_length x num_symbols matrix of p-values
    #         answer[4]   - a sequence_length x num_symbols matrix of frequency
    #                       differences
    def TSL.get_p_values(pos_data, neg_data, params)
        answer = Array.new

        # Get positive sequences
        pos        = get_seqs(pos_data, params['INPUT_KIND'])
        pos_seqs   = pos[0]
        pos_desc   = pos[1]
        pos_length = pos[2]
        error_code = pos[3]
        error_line = pos[4]

        # Check if positives are formatted properly
        if (1 == error_code)
            answer << 1 << error_line << nil << nil << nil
	    return answer
        end

        # Check non-uniform lengths
        if (2 == error_code)
            answer << 2 << error_line << nil << nil << nil
            return answer 
        end

        # Get negative sequences
        neg        = get_seqs(neg_data, params['INPUT_KIND'])
        neg_seqs   = neg[0]
        neg_desc   = neg[1]
        neg_length = neg[2]
        error_code = neg[3]
        error_line = neg[4]

        # Check if negatives are formatted properly
        if (1 == error_code)
            answer << 3 << error_line << nil << nil << nil
	    return answer
        end

        # Check non-uniform lengths
        if (2 == error_code)
            answer << 4 << error_line << nil << nil << nil
            return answer 
        end

        # Check that positives and negatives are the same length
        if (pos_length != neg_length)
            answer << 5 << nil << nil << nil << nil
            return answer
        end


        # Compute p-values
        p_values = Array.new
        frequencies = Array.new

        for i in 0...pos_length
            # Initialize symbol counters
            pos_counter = Hash.new
            neg_counter = Hash.new

            if ("A" == params['INPUT_KIND'])
                @@AMINOACIDS.each do |aa|
                    pos_counter[aa]=0
                    neg_counter[aa]=0
                end
            else
                @@NUCLEOTIDES.each do |na|
                    pos_counter[na]=0
                    neg_counter[na]=0
                end
            end

            # Count symbols
            pos_seqs.each do |seq|
                pos_counter[seq[i,1]] += 1  if (!space?(seq[i,1]))
            end
            neg_seqs.each do |seq|
                neg_counter[seq[i,1]] += 1  if (!space?(seq[i,1]))
            end

            # Compute frequency differences and prepare the p-value executable call
            frequencies << Array.new
            p_values << Array.new

            cmd ="../pvalue/pvalue "

            if ("binomial" == params['TEST'])
                cmd << "b"
            else
                cmd << "t"
            end


            pos_total = 0 
            neg_total = 0 

            if ("A" == params['INPUT_KIND'])
                @@AMINOACIDS.each do |aa|
                    if (pos_counter[aa]==pos_seqs.length && neg_counter[aa]==neg_seqs.length)
                        frequencies[i] << 1.0  # Conserved residue
                    else
                        frequencies[i] << (1.0 * pos_counter[aa]/pos_seqs.length) - (1.0 * neg_counter[aa]/neg_seqs.length)
                    end

                    cmd += " #{pos_counter[aa]} #{neg_counter[aa]}"

                    pos_total += pos_counter[aa]
                    neg_total += neg_counter[aa]
                    p_values[i] << 1.0
                end
            else
                @@NUCLEOTIDES.each do |na|
                    if (pos_counter[na]==pos_seqs.length && neg_counter[na]==neg_seqs.length)
                        frequencies[i] << 1.0  # Conserved residue
                    else
                        frequencies[i] << (1.0 * pos_counter[na]/pos_seqs.length) - (1.0 * neg_counter[na]/neg_seqs.length)
                    end

                    cmd += " #{pos_counter[na]} #{neg_counter[na]}"

                    pos_total += pos_counter[na]
                    neg_total += neg_counter[na]
                    p_values[i] << 1.0
                end 
            end

            if (pos_total>0 && neg_total>0)  # Symbols at this position, or only dashes/spaces?
                result = `#{cmd}`.split("\n")  # Call p-value executable

                for j in 0...result.length
                    if ("NAN" == result[j].upcase)
                        p_values[i][j] = *['7FF8000000000000'].pack('H*').unpack('G')  # One way to initialize NaN in Ruby
                    else
                        p_values[i][j] = result[j].to_f
                    end

                    if (params['BONFERRONI'])
                        if ("A" == params['INPUT_KIND'])
                            p_values[i][j] *= 20 * pos_length
                        else
                            p_values[i][j] *= 4  * pos_length
                        end
                    end
                end
            end
        end

        answer << 0 << nil << pos_length << p_values << frequencies
        return answer
    end


    # Builds arrays of heights.
    #
    # Input:  seq_length  - sequence length
    #         p_values    - a seq_length x num_symbols matrix of p-values
    #         frequencies - a seq_length x num_symbols matrix of frequency
    #                       differences
    #         params      - hash table containing parameters
    #
    # Output: An array containing the following:
    #         answer[0]  - conserved residues
    #         answer[1]  - enriched residues
    #         answer[2]  - depleted residues
    #         answer[3]  - maximum position height (for scaling)
    def TSL.get_heights(seq_length, p_values, frequencies, params)
        data_conserved = Array.new
        data_enriched  = Array.new
        data_depleted  = Array.new

        max_height = 0

        for i in 0...seq_length
            # Check for conserved symbols
            j = 0
            conserved = " "

            if ("A" == params['INPUT_KIND'])
                while (" "==conserved && j<@@AMINOACIDS.length)
                    conserved = @@AMINOACIDS[j]  if (1.0 == frequencies[i][j])
                    j += 1
                end
            else
                while (" "==conserved && j<@@NUCLEOTIDES.length)
                    conserved = @@NUCLEOTIDES[j]  if (1.0 == frequencies[i][j])
                    j += 1
                end
            end

            if (" " != conserved)
                data_conserved << conserved
                data_enriched << Array.new
                data_depleted << Array.new            
            else
                vert_enriched = Array.new
                vert_depleted = Array.new

                if ("A" == params['INPUT_KIND'])
                    for j in 0...(@@AMINOACIDS.length)
                        aa = @@AMINOACIDS[j]

                        if (p_values[i][j] < params['P_VALUE'])
                            height = frequencies[i][j]

                            if (params['FIXED_HEIGHT'])
                                if (height > 0)
                                    vert_enriched << "#{aa}0.05"
                                elsif (height < 0)
                                    vert_depleted << "#{aa}-0.05"
                                end
                            else
                                if (height > 0)
                                    vert_enriched << "#{aa}#{height}"
                                elsif (height < 0)
                                    vert_depleted << "#{aa}#{height}"
                                end
                            end
                        end
                    end
                else
                    for j in 0...(@@NUCLEOTIDES.length)
                        na = @@NUCLEOTIDES[j]

                        if (p_values[i][j] < params['P_VALUE'])
                            height = frequencies[i][j]

                            if (params['FIXED_HEIGHT'])
                                if (height > 0)
                                    vert_enriched << "#{na}0.05"
                                elsif (height < 0)
                                    vert_depleted << "#{na}-0.05"
                                end
                            else
                                if (height > 0)
                                    vert_enriched << "#{na}#{height}"
                                elsif (height < 0)
                                    vert_depleted << "#{na}#{height}"
                                end
                            end
                        end
                    end
                end
 
                # Sort enriched by height
                vert_enriched = vert_enriched.sort { |a, b|
                    a =~ /^(.{1})(\S+)$/
                    lettera = $1
                    heighta = $2
                    b =~ /^(.{1})(\S+)$/
                    letterb = $1
                    heightb = $2

                    # Compare by height first, then letter
                    if ((heighta <=> heightb) != 0)
                        heighta <=> heightb
                    else
                       letterb <=> lettera 
                    end
                }

                # Sort depleted by height
                vert_depleted = vert_depleted.sort { |a, b|
                    a =~ /^(.{1})(\S+)$/
                    lettera = $1
                    heighta = $2
                    b =~ /^(.{1})(\S+)$/
                    letterb = $1
                    heightb = $2

                    # Compare by height first, then letter
                    if ((heighta <=> heightb) != 0)
                        heighta <=> heightb
                    else
                        letterb <=> lettera
                    end
                }

                temp = 0
                vert_enriched.each do |x|
                    temp += x[1..x.length].to_f
                end
                max_height = temp  if (temp>max_height)


                temp = 0
                vert_depleted.each do |x|
                    temp -= x[1..x.length].to_f
                end
                max_height = temp  if (temp>max_height)

                data_conserved << " "
                data_enriched << vert_enriched
                data_depleted << vert_depleted
            end
        end

        answer = Array.new
        answer << data_conserved << data_enriched << data_depleted << max_height
        return answer
    end


    # Input:  input_data   -
    #         INPUT_KIND   - A for amino acids, N for nucleotides
    # 
    # Output: seqs         - array reference to sequence strings
    #         desc         - array reference to sequence names
    #         maxres_length - length of sequence
    #         goodlength   - 1 if all sequences have the same length, 0 else
    #         badline      - line number L where sequenceLength(L) != sequenceLength(other lines), 
    #                        else nil
    #         validformat  -
    #
    def TSL.get_seqs(input_data, input_kind)
        answer = Array.new

        if (0 == input_data.length)
            answer << nil << nil << nil << nil << nil << 1
            return answer
        end  
 
        # skip comments and empty lines
        while (input_data[0] =~ /^\s*\#/ || input_data[0] =~ /^\s*$/)
            input_data.shift
        end

        return get_seqs_FASTA(input_data, input_kind)   if format_FASTA?(input_data)
        return get_seqs_CLUSTAL(input_data, input_kind) if format_CLUSTAL?(input_data)
	return get_seqs_FLAT(input_data, input_kind)    if format_FLAT?(input_data)

        answer << nil << nil << nil << nil << nil << 1
        return answer
    end


    #
    # Returns true if the argument is in FASTA format
    #
    def TSL.format_FASTA?(input_data)
        return true  if (">" == input_data[0][0,1])
        return false
    end


    #
    # Returns true if the argument is in CLUSTAL format
    #
    def TSL.format_CLUSTAL?(input_data)
        # if it looks like CLUSTAL W (version) ... , then it must be clustal
        return true  if (input_data[0] =~ /^\s*CLUSTAL/)

        # CLUSTAL looks like: "name        seq"
        return true  if (input_data[0] =~ /^\s*(\S+)\s+(\S+)\s*$/)
 
        return false
    end


    #
    # Returns true if the argument is a flat file. 
    #
    def TSL.format_FLAT?(input_data)
        return true  if (input_data[0] =~ /^[a-zA-Z\-]+\s*$/)
        return false
    end


    #
    # Masks all but ACTGU for nucleotides, and 
    # ACDEFGHIKLMNPQRSTVWY for amino acids. 
    #
    def TSL.filter_char(ch, input_kind)
        char = ""

        if ("A" == input_kind)
            char = @@AMINOACIDS.include?(ch) ? ch : "-"
        else
            char = @@NUCLEOTIDES.include?(ch) ? ch : "-"
        end

        return char
    end


    #
    # Get sequences from a flat sequence file.
    #
    # Output: seqs         - array of sequences
    #         desc         - nil (it's a flat file)
    #         length       - sequence length
    #         error_code   -
    #         error_line   -
    #
    def TSL.get_seqs_FLAT(input_data, input_kind)
        answer       = Array.new
        sequences    = Array.new
        counter      = 0
        first_length = nil

        input_data.each do |input| 
            input = input.sub(/\s+$/, "")

            chars = input.split(//)

            chars.each do |ch|
                sequences[counter]  = "" if sequences[counter]==nil
                sequences[counter] += filter_char(ch, input_kind)
            end

            if counter==0
                first_length = input.length
            elsif first_length != input.length             # different number of residues, so complain
                answer << nil << nil << nil << 2 << input  # 0 for not same length, input is wrong
                return answer
            end
	
            counter += 1
        end

        answer << sequences << Array.new << first_length << 0 << nil
        return answer
    end


    #
    # Get sequences from a file in CLUSTAL format.
    #
    def TSL.get_seqs_CLUSTAL(input_data, input_kind)
        answer     = Array.new
        sequences  = Array.new
        desc       = Array.new
        seq_count  = 0
        res_length = 0
        name       = ""
        seq        = ""

        prevline_length = 0
        line_length     = 0

        input_data.each do |input|
            input = input.sub(/\s+$/, "")

            # skip if it is a comment character -- first character is "#"
            next if input =~ /^\s*\#/

            # skip if it is a CLUSTAL W header line
            next if input =~ /^\s*CLUSTAL/

            # if spaces or just "*" and "." and ":"
            if input =~ /^[\*\.\:\s]*$/ 
                seq_count = 0
                prevline_length = 0
                next
            end	

            input =~ /^\s*(\S+)\s+(\S+)\s*$/
            name = $1
            seq  = $2

            # add new entry
            if (desc[seq_count] == nil)
                desc[seq_count]      = name
                sequences[seq_count] = ""
            end

            chars = seq.split(//)

            chars.each do |ch|
                # all sequences have same residue length, so only count first one
                res_length +=1 if (seq_count == 0)
                line_length += 1
	    
                sequences[seq_count] += filter_char(ch, input_kind)
            end

            if seq_count == 0
                prevline_length = line_length
            elsif prevline_length != line_length  # different number of residues, so complain
                answer << nil << nil << nil << 2 << name
                return answer
            end

            line_length =  0
            seq_count  += 1
        end 

        answer << sequences << desc << res_length << 0 << nil
        return answer
    end


    #
    # Get sequences from a file in FastA format.
    # 
    def TSL.get_seqs_FASTA(input_data, input_kind)
        answer        = Array.new
        sequences     = Array.new
        desc          = Array.new
        count         = -1
        newElem       = false
        res_length    = 0
        maxres_length = 0
    
        goodlength = ""
        currline   = ""
        prevline   = ""

        input_data.each do |input| 

            input = input.sub(/\s+$/, "")

            # skip if it is a comment character -- first character is "#"
            next if input =~ /^\s*\#/

            # skip all lines that are all spaces
            next if input =~ /^\s*$/

            input = input.strip 

            if input =~ />/
                currline = input 
                desc[desc.length] = $1 if input =~ />\s*(.+)$/

                if !newElem 		
                    count   += 1
                    newElem = true
                end
            else
	        if newElem
		    maxres_length = res_length if maxres_length == 0
                    if maxres_length != 0 && maxres_length != res_length
                        answer << nil << nil << nil << 2 << prevline
                        return answer
                    end

                    maxres_length = res_length
                    res_length = 0
                end 

                chars = input.split(//)
 
                chars.each do |ch|
                    res_length += 1

                    sequences[count] = "" if newElem
                    sequences[count] += filter_char(ch, input_kind)

                    newElem = false 
                end
                prevline = currline if currline =~ />/
            end
        end

        # check if last is biggest
        if maxres_length != 0 && maxres_length != res_length
            answer << nil << nil << nil << 2 << prevline
            return answer
        end

        answer << sequences << desc << maxres_length << 0 << nil
        return answer
    end


    #
    # Checks if the argument is space or a gap.
    #
    def TSL.space?(arg)
        return arg =~ /[ \-]/
    end




end
