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
