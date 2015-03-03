function rev(string, x,i) {
		for(i=length;i!=0;i--)
				x=(x substr(string,i,1))
		return x
}

function abs(value) {
    return (value<0?-value:value);
}

function div(denom, divd) {
    return (divd == 0? "NA" : denom/divd);
}

function dnarev(string,        chars, dic, new_string) {
    dic["A"] = "T";
    dic["T"] = "A";
    dic["G"] = "C";
    dic["C"] = "G";

    new_string = "";
    split(string, chars, "");
    for (i=length(chars);i>=0;i--) {
        new_string = new_string dic[chars[i]];
    }

    return new_string;
}

# join fields by specifying a range
# for example, from $3 to $10, one can use joinr(3,10)
function joinr(start, end,     result, i) {

    result = ""
    for (i=start; i<=end; ++i) {
	if (i != start) {
	    result = result "\t";
	}
	result = result $i;
    }
    return result
}

function join(array, start, end, sep,    result, i)
{
    if (sep == "")
       sep = " "
    else if (sep == SUBSEP) # magic value
       sep = ""
    result = array[start]
    for (i = start + 1; i <= end; i++)
        result = result sep array[i]
    return result
}

# escape the string when the string is
# used as a pattern in regular expression
# i.e., \ needs to be replaced by \\
function escape(pat,   safe) {
  safe = pat
  gsub(/[][^$.*?+{}\\()|]/, "\\\\&", safe)
  return safe
}