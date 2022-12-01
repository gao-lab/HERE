# Util.rb - Utility class for processing command line options.
#
# Vladimir Vacic <vladimir@cs.ucr.edu>
# Algorithms and Computational Biology Lab 
# University of California, Riverside
#
# Aug-18-2005


class Util

    #
    # Usage: opts = getopts("a:bc", ARGV)   -a takes arg. -b & -c not
    #
    def Util.getopts(switches)
        opts = Hash.new

        while (!ARGV.empty? && ARGV[0] =~ /^-(.)(.*)/)
            first = $1 
            rest  = $2

            pos = switches.index(first)

            if pos != nil
                if switches[pos+1, 1] == ":"
                    ARGV.shift
                    rest = ARGV.shift if rest == "" 

                    if opts[first] == nil
                        opts[first] = rest
                    else
                        opts[first] += " #{rest}"
                    end 
                else
                    opts[first] = true
		    if rest == ""
                        ARGV.shift
		    else
                        ARGV[0] = "-#{rest}"
		    end
                end	
            else
                puts "Unknown option: #{first}\n"
                if rest != ""
                    ARGV[0] = "-#{rest}"
                else
                    ARGV.shift
                end
            end
        end

        return opts
    end
end
