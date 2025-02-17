import os
import glob
import time

def mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)
    return os.path.abspath(path)
def rmfile(file):
    if os.path.exists(file):
        os.remove(file)

def backup_std_err(sample, sub_log_dir):
    # 备份日志文件
    sub_log = f"{sub_log_dir}/{sample}.run.log"
    flag = 1
    if not os.path.exists(sub_log):
        logf=open(sub_log,'w')
    else:
        cmd_list = glob.glob(f"{sub_log_dir}/{sample}.run.std.bak*")
        num = flag + len(cmd_list)
        bak_log = f"{sub_log_dir}/{sample}.run.log.bak{num}"
        os.rename(sub_log, bak_log)
        logf = open(sub_log, 'w')
    return logf
def startTime():
    start_time = time.time()
    print(time.strftime('Job starts at %Y/%m/%d %H:%M:%S', time.localtime(start_time)))
    return start_time


def spendTime(start_time):
    end_time = time.time()
    print(time.strftime('Job ends at %Y/%m/%d %H:%M:%S', time.localtime(end_time)))
    spend = end_time - start_time
    print(f"Real time: {spend:.4f}")
    h = int(spend / 3600)
    remain = spend % 3600
    m = int(remain / 60)
    remain = remain % 60
    s = int(remain)
    return f"{h}h{m}m{s}s"


