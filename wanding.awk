function dnarev(string) {
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
# used as a pattern
function escape(pat,   safe) {
  safe = pat
  gsub(/[][^$.*?+{}\\()|]/, "\\\\&", safe)
  return safe
}
