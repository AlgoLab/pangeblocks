import sys 

# return ['user_time', 'sys_time', 'max_mem', 'wall_clock'] from log /usr/bin/time

def main(path_log):
    
    # load log
    with open(path_log) as fp: 
        res=dict()
        for line in fp.readlines():
            line=line.strip()
            tokens = line.split(sep=":")

            if tokens[0] == "User time (seconds)":
                res['user_time'] = float(tokens[1].lstrip())
            if tokens[0] == "System time (seconds)":
                res['sys_time'] = float(tokens[1].lstrip())
            if tokens[0] == "Maximum resident set size (kbytes)":
                res['max_mem'] = int(tokens[1].lstrip())
            if tokens[0] == "Elapsed (wall clock) time (h":
                tot = 0.0
                for x in tokens[4:]:
                    tot = tot*60 + float(x.lstrip())
                res['wall_clock'] = tot

    return res

if __name__ == "__main__":

    path_log=sys.argv[1]
    res=main(path_log)  
    print(res)