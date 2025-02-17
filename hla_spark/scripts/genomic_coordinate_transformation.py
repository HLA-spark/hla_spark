from collections import OrderedDict
import numpy as np
import pandas as pd

def retrive_ref(df_pos, chr,gene):
    ref = df_pos.loc[gene].values
    if chr in df_pos.index.values:
        keys = df_pos.loc[chr].values
    else:
        keys = ref
    df_ref = pd.DataFrame(columns=range(0, len(keys)))
    df_ref.loc['ref', :] = ref
    df_ref.loc['ref_before', :] = keys
    df_ref = df_ref.transpose()
    df_ref['ref_base'] = df_ref.apply(lambda x: str(x['ref']).split('.')[0], axis=1)
    df_ref['ref_pos'] = df_ref.apply(lambda x: str(x['ref']).split('.')[1], axis=1)
    return df_ref


def indel_judge(tokens):
    info=tokens[7]
    if 'INDEL' in info:
        return 'INDEL'
    else:
        return 'SNP'

def s01_line2info(tokens):
    l_pos = tokens[1]
    l_ref = tokens[3]
    l_alt=tokens[4].split(',')
    l_info=tokens[7]
    l_dp=int(l_info.split('DP=')[1].split(';')[0])
    l_i16=l_info.split('I16=')[1].split(';')[0].split(',')
    for i in range(0,len(l_i16)):
        l_i16[i]=str(int(eval(l_i16[i])))

    l_qs=l_info.split('QS=')[1].split(';')[0].split(',')

    formt = tokens[9].split(':')

    l_dp_qc = int(formt[1])
    l_ad = formt[2].split(',')
    l_pl = formt[0].split(',')

    if '<*>' not in l_alt:
        ori_n=len(l_ad)
        new_n=ori_n+1
        ori_n_pl=int(0.5*ori_n*(ori_n+1))
        new_n_pl=int(0.5*new_n*(new_n+1))
        l_alt.append('<*>')
        l_qs.append('0')
        l_ad.append('0')
        l_pl.extend((new_n_pl-ori_n_pl)*['255'])

    return l_pos, l_ref, l_alt, l_dp,l_i16,l_qs,l_pl,l_dp_qc, l_ad

def s02_info_conversion(l_pos, l_ref, l_alt, l_dp,l_i16,l_qs,l_pl,l_dp_qc, l_ad, df_ref):
    r_key = l_ref + '.' + l_pos
    ref_base, pos = s03_get_real_refbase_pos(r_key, df_ref)
    if l_ref==ref_base:
        return pos, ref_base, ','.join(l_alt), l_dp,','.join(l_i16),','.join(l_qs), ','.join(l_pl), l_dp_qc,','.join(l_ad)
    else:
        if ref_base in l_alt:# ALT中有与ref_base相同的碱基突变
            same_index=s04_get_same_index(ref_base,l_alt)#same_index为l_alt中与ref_base相同的碱基的坐标
            alt_base_list,i16,qs_list,pl_list,ad_list=s05_ref_in_alt(l_ref,l_alt,same_index,l_i16,l_qs,l_pl,l_ad)

        else:#alt中没有与ref_base相同的碱基突变，ref_base作为新的ref，之前的ref加入首个alt中
            alt_base_list, i16, qs_list, pl_list, ad_list = s09_ref_not_in_alt(l_ref, l_alt, l_i16, l_qs, l_pl,l_ad,pos)

        order_alt_base_list, order_qs_list, order_pl_list, order_ad_list = s13_sort_alt_ad(alt_base_list, qs_list, pl_list,ad_list)

        return pos, ref_base, order_alt_base_list, l_dp, i16, order_qs_list, order_pl_list, l_dp_qc, order_ad_list

def s03_get_real_refbase_pos(ref_before, df_ref):  # 通过转变前ref位置+碱基获取对应的转变后ref位置+碱基
    df_relation = df_ref[df_ref['ref_before'] == ref_before]
    real_refs = df_relation['ref'].values[0].split('.')
    real_base = real_refs[0]
    real_pos = real_refs[1]

    if real_base == '-':
        df_relation = df_ref[(df_ref['ref_pos'] == real_pos) & (df_ref['ref_base'] != '-')]
        if len(df_relation['ref'].values)<1:
            return real_base, real_pos
        real_refs = df_relation['ref'].values[0].split('.')
        real_ref_before = df_relation['ref_before'].values[0].split('.')
        inc = int(ref_before.split('.')[1]) - int(real_ref_before[1])
        real_pos = real_refs[1] + '_' + str(inc)

    return real_base, real_pos

def s04_get_same_index(ref_base,l_alt):
    for i in range(0,len(l_alt)):
        if ref_base==l_alt[i]:
            return i

def s05_ref_in_alt(l_ref,l_alt,same_index,l_i16,l_qs,l_pl,l_ad):
    alt_base_list=[l_ref]
    l_alt.pop(same_index)
    alt_base_list.extend(l_alt)
    l_ad_int=[]
    for i in l_ad:
        l_ad_int.append(int(i))
    #l_ad最后变动
    if len(l_ad)==2:#ref_alt交换位置
        i16=s06_i16_exchange_single(l_i16)
    else:
        ratio=l_ad_int[1+same_index]/(sum(l_ad_int[1:]))
        i16=s07_i16_exchange_multi(l_i16,ratio)

    qs_list=[l_qs[same_index+1]]
    l_qs.pop(same_index+1)
    qs_list.extend(l_qs)

    pl_list=s08_pl_exchange(l_pl,same_index,len(l_ad))

    ad_list = [l_ad[1 + same_index]]
    l_ad.pop(same_index + 1)
    ad_list.extend(l_ad)
    return alt_base_list,i16,qs_list,pl_list,ad_list

def s06_i16_exchange_single(l_i16):
    index=[2,3,0,1,6,7,4,5,10,11,8,9,14,15,12,13]
    i16=[]
    for i in index:
        i16.append(l_i16[i])
    return ','.join(i16)

def s07_i16_exchange_multi(l_i16,p):
    index1 = [2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13]
    index2=list(range(0,16))
    ratio1=[p,p,1,1,p,p,1,1,p,p,1,1,p,p,1,1]
    ratio2=[0,0,1-p,1-p,0,0,1-p,1-p,0,0,1-p,1-p,0,0,1-p,1-p]
    i16=[]
    i16_int=[]
    for i in l_i16:
        i16_int.append(int(i))
    for i in index2:
        i16.append(str(int((ratio1[i]*i16_int[index1[i]]+ratio2[i]*i16_int[index2[i]]))))
    return ','.join(i16)

def s08_pl_exchange(l_pl,same_index,allelic_num):
    same_index=same_index+1
    n=allelic_num
    n_pl=int(0.5*n*(n+1))
    pl_type_order = ['0/0', '0/1', '1/1', '0/2','1/2', '2/2', '0/3', '1/3','2/3', '3/3', '0/4','1/4', '2/4','3/4','4/4','0/5','1/5','2/5','3/5','4/5','5/5']
    pl_type=pl_type_order[0:n_pl]
    l_pl_dict=dict()
    for i in range(0,len(pl_type)):
        l_pl_dict[pl_type[i]]=l_pl[i]
    index_change_dict=dict()
    index_st=2
    for i in range(0,n):
        if i==0:
            index_change_dict[i]=1
        elif i==same_index:
            index_change_dict[i]=0
        else:
            index_change_dict[i]=index_st
            index_st+=1

    unsort_pl_dict =dict()
    for key_before in l_pl_dict.keys():
        tokens = key_before.split('/')
        after1 = index_change_dict.get(int(tokens[0]))
        after2 = index_change_dict.get(int(tokens[1]))
        new_key = '/'.join([str(min(after1, after2)), str(max(after1, after2))])
        unsort_pl_dict[new_key]=l_pl_dict.get(key_before)
    new_pl_list=[]
    for i in pl_type:
        new_pl_list.append(unsort_pl_dict.get(i))
    return new_pl_list

def s09_ref_not_in_alt(l_ref, l_alt, l_i16, l_qs, l_pl,l_ad,pos):
    # alt中没有与ref_base相同的碱基突变，ref_base作为新的ref，之前的ref加入首个alt中
    alt_base_list=[l_ref]
    alt_base_list.extend(l_alt)

    i16=s10_i16_exchange_newref(l_i16)

    qs_list=['0']
    qs_list.extend(l_qs)

    ad_list=['0']
    ad_list.extend(l_ad)
    n=len(ad_list)
    pl_list=s11_pl_exchange_newref(l_pl,n-1,n,pos)

    return alt_base_list,i16,qs_list,pl_list,ad_list

def s10_i16_exchange_newref(l_i16):
    ori_ref_index=[0,1,4,5,8,9,12,13]
    ori_non_ref_index=[2,3,6,7,10,11,14,15]
    i16_int=[]
    flag=0
    j=0
    for i in range(0,8):
        if flag==0:
            i16_int.extend([0,0])
            flag=1
        else:
            i16_int.extend([int(l_i16[ori_ref_index[j]])+int(l_i16[ori_non_ref_index[j]]),int(l_i16[ori_ref_index[j+1]])+int(l_i16[ori_non_ref_index[j+1]])])
            j=j+2
            flag=0
    i16=[]
    for i in i16_int:
        i16.append(str(i))
    return ','.join(i16)

def s11_pl_exchange_newref(l_pl,ori_n,new_n,pos):
    ori_n_pl=int(0.5*ori_n*(ori_n+1))
    new_n_pl=int(0.5*new_n*(new_n+1))
    pl_type_order = ['0/0', '0/1', '1/1', '0/2','1/2', '2/2', '0/3', '1/3','2/3', '3/3', '0/4','1/4', '2/4','3/4','4/4','0/5','1/5','2/5','3/5','4/5','5/5']
    ori_pl_type=pl_type_order[0:ori_n_pl]
    new_pl_type=pl_type_order[0:new_n_pl]
    ori_pl_dict=dict()
    for i in range(0,ori_n_pl):
        ori_pl_dict[ori_pl_type[i]]=l_pl[i]

    change_dict=OrderedDict()
    new_known_index_set=set()
    for key_before in ori_pl_type:
        tokens=key_before.split('/')
        after1=int(tokens[0])+1
        after2=int(tokens[1])+1
        new_key='/'.join([str(min(after1,after2)),str(max(after1,after2))])
        change_dict[new_key]=ori_pl_dict.get(key_before)
        new_known_index_set.add(str(after1))
        new_known_index_set.add(str(after2))
    change_dict=s12_heter_new_pl_fill(change_dict,new_n,new_pl_type,new_known_index_set,pos)
    new_pl_list=[]
    for i in new_pl_type:
        new_pl_list.append(change_dict.get(i))

    return new_pl_list

def s12_heter_new_pl_fill(change_dict,new_n,new_pl_type,know_index_set,pos):
    null_index = new_n - 1
    for i in new_pl_type:
        if i not in change_dict.keys():
            tokens = i.split('/')
            index1 = int(tokens[0])
            index2 = int(tokens[1])
            if str(index1) not in know_index_set:
                index1 = null_index
            if str(index2) not in know_index_set:
                index2 = null_index
            rep_key = '/'.join([str(min(index1, index2)), str(max(index1, index2))])
            change_dict[i] = change_dict.get(rep_key)
    return change_dict

def s13_sort_alt_ad(alt_base_list, qs_list, pl_list,ad_list):
    ref_depth=ad_list[0]
    ref_qs=str(qs_list[0])
    dict_base_ad = dict()
    for i in range(0, len(alt_base_list)):
        dict_base_ad[alt_base_list[i]] = ad_list[i+1]

    dict_sort = dict(sorted(dict_base_ad.items(), key=lambda x: x[1], reverse=True))

    new_altbase_list = list(dict_sort.keys())
    new_altbase = ','.join(new_altbase_list)

    new_ad_list = [ref_depth]
    new_ad_list.extend(dict_sort.values())
    new_ad = ','.join(new_ad_list)

    ori_altbase_index_dict=dict()
    for i in range(0,len(alt_base_list)):
        ori_altbase_index_dict[alt_base_list[i]]=i+1

    new_altbase_index_dict=dict()
    for i in range(0,len(new_altbase_list)):
        new_altbase_index_dict[new_altbase_list[i]]=i+1


    ##qs
    new_qs_list=[ref_qs]
    new_qs_list.extend(len(new_altbase_list)*['0'])
    for i in range(1,len(qs_list)):
        altbase=alt_base_list[i-1]
        qs=qs_list[i]
        new_index=new_altbase_index_dict.get(altbase)
        new_qs_list[new_index]=qs
    new_qs = ','.join(new_qs_list)

    ##pl
    n=len(ad_list)
    n_pl=int(0.5*n*(n+1))
    pl_type_order = ['0/0', '0/1', '1/1', '0/2','1/2', '2/2', '0/3', '1/3','2/3', '3/3', '0/4','1/4', '2/4','3/4','4/4','0/5','1/5','2/5','3/5','4/5','5/5']
    pl_type=pl_type_order[0:n_pl]
    convert_ori_pl_dict=dict()
    for i in range(0,len(pl_type)):
        new1=0
        new2=0
        tokens=pl_type[i].split('/')
        before1=tokens[0]
        before2=tokens[1]
        if before1!='0':
            base1=alt_base_list[int(before1)-1]
            new1=new_altbase_index_dict.get(base1)
        if before2!='0':
            base2=alt_base_list[int(before2)-1]
            new2=new_altbase_index_dict.get(base2)
        news=[str(min(new1,new2)),str(max(new1,new2))]
        convert_ori_pl_dict['/'.join(news)]=pl_list[i]
    new_pl_list=[]
    for i in pl_type:
        new_pl_list.append(convert_ori_pl_dict.get(i))
    new_pl = ','.join(new_pl_list)
    return new_altbase, new_qs, new_pl, new_ad

def s14_heter(ref_base, add_altbase, add_dp,add_i16,add_qs,add_pl,add_dp_qc, add_ad, values,pos):#values:['REF', 'ALT', 'DP','I16','QS','PL','DP_qc', 'AD']

    ori_altbase_list = values[1].split(',')
    ori_dp_int = int(values[2])
    ori_i16_list=values[3].split(',')
    ori_qs_list=values[4].split(',')
    ori_dp_qc_int=int(values[6])
    ori_ad_list = values[7].split(',')

    new_dp=str(ori_dp_int+int(add_dp))
    new_dp_qc=str(ori_dp_qc_int+int(add_dp_qc))

    new_i16,baseQ=s15_i16_add(ori_i16_list,add_i16.split(','))

    add_ad_list = add_ad.split(',')
    add_altbase_list = add_altbase.split(',')
    add_qs_list=add_qs.split(',')

    new_ref_depth = int(ori_ad_list[0]) + int(add_ad_list[0])

    new_alt_depth_dict = dict()  #记录最终alt的深度信息
    ori_alt_base_index_dict=dict() #记录原始alt的位置信息

    for i in range(1, len(ori_ad_list)):
        ori_alt_base_index_dict[ori_altbase_list[i - 1]] = i
        new_alt_depth_dict[ori_altbase_list[i - 1]] = int(ori_ad_list[i])


    add_alt_base_index_dict=dict() #记录新增alt的位置信息
    for i in range(1, len(add_ad_list)):
        add_alt_base_index_dict[add_altbase_list[i - 1]] = i
        if add_altbase_list[i - 1] in new_alt_depth_dict.keys():
            od = new_alt_depth_dict.get(add_altbase_list[i - 1])
            new_alt_depth_dict[add_altbase_list[i - 1]] = int(od) + int(add_ad_list[i])
        else:
            new_alt_depth_dict[add_altbase_list[i-1]]=int(add_ad_list[i])

    ##根据深度对alt排序（ref单独计算）
    new_alt_depth_dict['<*>']=-1

    dict_sort = dict(sorted(new_alt_depth_dict.items(), key=lambda x: x[1], reverse=True))
    dict_sort['<*>']=0
    new_altbase_list=list(dict_sort.keys())#坐标+1为对应碱基的位置信息
    new_alt_depth_list=[]
    for i in dict_sort.values():
        new_alt_depth_list.append(str(i))
    new_altbase = ','.join(new_altbase_list)
    new_ad_list=[str(new_ref_depth)]
    new_ad_list.extend(new_alt_depth_list)
    new_ad=','.join(new_ad_list)

    new_altbase_index_dict=dict()
    for i in range(0,len(new_altbase_list)):
        new_altbase_index_dict[new_altbase_list[i]]=i+1
    ##原始指标根据排序后alt进行位置调整和基因型替换 qs和pl
    ##新增指标根据排序后alt进行位置调整和基因型替换 qs和pl
    ##原始和新增指标合并
    ori_dp_ratio=0
    add_dp_ratio=0
    if int(new_dp_qc)>0:
        ori_dp_ratio = 1.0 * ori_dp_qc_int / int(new_dp_qc)
        add_dp_ratio = 1.0 * int(add_dp_qc) / int(new_dp_qc)

    ##qs
    convert_ori_qs_list=s18_heter_qs_convert(ori_altbase_list,new_altbase_index_dict,ori_qs_list)#float
    convert_add_qs_list=s18_heter_qs_convert(add_altbase_list,new_altbase_index_dict,add_qs_list)#float

    new_qs_list = []
    for i in range(0, len(convert_ori_qs_list)):
        x = convert_ori_qs_list[i]
        y = convert_add_qs_list[i]
        new_qs_list.append(ori_dp_ratio * x + add_dp_ratio * y)

    new_qs_format_list=[]
    qs_list=[] #float型format qs list
    qs_sum=sum(new_qs_list)
    if qs_sum>0:
        for j in new_qs_list:
            new_qs_format_list.append(str(1.0*round(j/qs_sum,6)))
            qs_list.append(1.0*round(j/qs_sum,6))
    else:
        for j in new_qs_list:
            new_qs_format_list.append(str(j))
            qs_list.append(j)
    new_qs=','.join(new_qs_format_list)

    ad_list=[]
    for i in new_ad_list:
        ad_list.append(int(i))

    new_pl=s16_heter_pl(baseQ,qs_list,ad_list,pos)

    return [ref_base, new_altbase, new_dp,new_i16,new_qs,new_pl,new_dp_qc, new_ad]

def s15_i16_add(ori_list,add_list):
    new_list=[]
    for i in range(0,len(ori_list)):
        new_list.append(str(int(ori_list[i])+int(add_list[i])))
    return ','.join(new_list),int(new_list[4])+int(new_list[6])

def s16_heter_pl(baseQ,qs_list,ad_list,pos):
    ##pl
    new_n = len(ad_list)
    new_n_pl = int(0.5 * new_n * (new_n + 1))
    depth = sum(ad_list)
    if depth == 0:
        return ','.join(new_n_pl* ['0'])

    pl_type_order = ['0/0', '0/1', '1/1', '0/2', '1/2', '2/2', '0/3', '1/3', '2/3', '3/3', '0/4', '1/4', '2/4', '3/4',
                     '4/4', '0/5', '1/5', '2/5', '3/5', '4/5', '5/5']
    new_pl_type = pl_type_order[0:new_n_pl]

    # pl合并
    e_list=[]
    p_list=[]
    for index in range(0,new_n):
        i=qs_list[index]
        j=ad_list[index]
        if j==0:
            E=0
        else:
            Q=1.0*baseQ*i/j
            E=(0.1)**(Q/10)
        e_list.append(E)
        if E==0:
            p_list.append(0)
        else:
            p_list.append(1-E)
    gl_list=[]
    for i in new_pl_type:
        tokens=i.split('/')
        h0=int(tokens[0])
        h1=int(tokens[1])
        gl_i=1
        for h in range(0,new_n):
            ad_h=ad_list[h]
            if ad_h==0:
                continue
            if h0 == h1:
                if h==h0:
                    gl_i *= p_list[h] ** ad_h
                else:
                    gl_i*=e_list[h]**ad_h
            else:
                if h==h0:
                    gl_i*=(0.5*p_list[h]+0.5*e_list[h1])**ad_h
                elif h==h1:
                    gl_i*=(0.5*p_list[h]+0.5*e_list[h0])**ad_h
                else:
                    gl_i*=(0.5*e_list[h1]+0.5*e_list[h0])**ad_h

        gl_list.append(gl_i)
    pl_list = []
    min_gl=1
    for gl_i in gl_list:
        if (gl_i>0) and (min_gl>gl_i):
            min_gl=gl_i

    for gl_i in gl_list:
        if gl_i==0:
            pl_i=-10*np.log10(min_gl)
            pl_i=pl_i+255
            if pl_i<255:
                pl_i=255
        else:
            pl_i=-10*np.log10(gl_i)
        pl_list.append(pl_i)
    norm_pl_list1=[]
    for i in pl_list:
        norm_pl_list1.append(i-min(pl_list))
    norm_coef = 1.0 * max(norm_pl_list1) / 255
    norm_pl_list=[]
    for i in norm_pl_list1:
        norm_pl_list.append(str(int(round(1.0*i/norm_coef))))
    new_pl=','.join(norm_pl_list)
    return new_pl

def s17_heter_target_pl_convert(target_pl_type,target_pl_list,target_altbase_list,new_altbase_index_dict):
    convert_target_pl_dict = dict()
    target_know_index_set = set()
    for i in range(0, len(target_pl_type)):
        tokens = target_pl_type[i].split('/')
        val = int(target_pl_list[i])
        new_index1 = 0
        new_index2 = 0
        if tokens[0] != '0':
            before_base1 = target_altbase_list[int(tokens[0]) - 1]
            new_index1 = new_altbase_index_dict.get(before_base1)
        if tokens[1] != '0':
            before_base2 = target_altbase_list[int(tokens[1]) - 1]
            new_index2 = new_altbase_index_dict.get(before_base2)
        target_know_index_set.add(str(new_index1))
        target_know_index_set.add(str(new_index2))
        new_genotype = '/'.join([str(min(new_index1, new_index2)), str(max(new_index1, new_index2))])
        convert_target_pl_dict[new_genotype] = val
    return convert_target_pl_dict,target_know_index_set

def s18_heter_qs_convert(target_altbase_list,new_altbase_index_dict,target_qs_list):#reture float list
    n=len(new_altbase_index_dict.keys())
    convert_target_qs_list=[float(target_qs_list[0])]
    convert_target_qs_list.extend(n*[0])
    for i in range(0,len(target_altbase_list)):
        alt=target_altbase_list[i]
        new_index=new_altbase_index_dict.get(alt)
        convert_target_qs_list[new_index]=float(target_qs_list[i+1])
    return convert_target_qs_list

def s19_heter_pl_combine(convert_ori_pl_dict,convert_add_pl_dict,new_pl_type,new_n,pos):
    n_pl=len(new_pl_type)
    convert_ori_pl_list=list(convert_ori_pl_dict.values())
    convert_add_pl_list=list(convert_add_pl_dict.values())

    if sum(convert_ori_pl_list)==0 and sum(convert_add_pl_list)==0:
        new_pl=n_pl*['0']
        return ','.join(new_pl)
    elif sum(convert_ori_pl_list)==0:
        new_pl=[]
        for i in convert_add_pl_list:
            new_pl.append(str(i))
        return ','.join(new_pl)
    elif sum(convert_add_pl_list)==0:
        new_pl=[]
        for i in convert_ori_pl_list:
            new_pl.append(str(i))
        return ','.join(new_pl)
    else:
        ori_pl_inc=255-max(convert_ori_pl_list)
        add_pl_inc=255-max(convert_add_pl_list)
        ori_gl_dict=dict()
        add_gl_dict=dict()
        for index in range(0,n_pl):
            i=convert_ori_pl_list[index]
            new_i=ori_pl_inc+i
            i_gl=10**(-new_i/10)
            ori_gl_dict[new_pl_type[index]]=i_gl
        for index in range(0,n_pl):
            i=convert_add_pl_list[index]
            new_i=add_pl_inc+i
            i_gl=10**(-new_i/10)
            add_gl_dict[new_pl_type[index]]=i_gl

        #new_alleles=list(range(0,new_n))
        ori_allele_p_list=new_n*[0]
        add_allele_p_list=new_n*[0]
        for i in ori_gl_dict.keys():
            val=ori_gl_dict.get(i)
            tokens=i.split('/')
            for j in tokens:
                j_int=int(j)
                ori_allele_p_list[j_int]+=val/2
        for i in add_gl_dict.keys():
            val=add_gl_dict.get(i)
            tokens=i.split('/')
            for j in tokens:
                j_int=int(j)
                add_allele_p_list[j_int]+=val/2
        allele_p_list=list(np.add(ori_allele_p_list,add_allele_p_list))
        min_allele_p=min(allele_p_list) #每一项减去最小值，使最小值概率为0
        for i in range(0,new_n):
            allele_p_list[i]-=min_allele_p
        ##进行标准化
        sum_allele_p=sum(allele_p_list)
        for i in range(0,new_n):
            allele_p_list[i]=1.0*allele_p_list[i]/sum_allele_p

        new_gl_list=n_pl*[0]
        for i in range(0,n_pl):
            i_type=new_pl_type[i]
            tokens=i_type.split('/')
            index1=int(tokens[0])
            index2=int(tokens[1])
            if index1==index2:
                new_gl_list[i]=allele_p_list[index1]*allele_p_list[index2]
            else:
                new_gl_list[i]=2*allele_p_list[index1]*allele_p_list[index2]

        new_pl_list=n_pl*[255]
        for i in range(0,n_pl):
            i_gl=new_gl_list[i]
            if i_gl!=0:
                new_pl_list[i]=-10*np.log10(i_gl)

        min_pl=min(new_pl_list)
        norm_pl_list=[]
        for i in range(0,n_pl):
            norm_pl_list.append(str(int(round(new_pl_list[i]-min_pl))))
        return ','.join(norm_pl_list)

def mpileupConvertToREF(gene_pos_matrix, mpile_vcf, ref_vcf,chrom,sample):
    df_pos = pd.read_csv(gene_pos_matrix, sep='\t', index_col=0, dtype=str)
    gene=df_pos.index[0]
    chr = ''
    df_ref_xls = pd.DataFrame(columns=['REF', 'ALT', 'DP','I16','QS','PL','DP_qc', 'AD'])
    vcf = open(mpile_vcf, 'r')
    lines = vcf.readlines()
    for line in lines:
        line = line.strip()
        if '#' in line:
            continue
        else:
            tokens = line.split('\t')
            var_type=indel_judge(tokens)
            if var_type=='INDEL': ##仅保留单个碱基的位置对应关系
                continue
            if tokens[0] != chr:
                chr = tokens[0]
                df_ref = retrive_ref(df_pos, chr,gene)
            l_pos, l_ref, l_alt, l_dp,l_i16,l_qs,l_pl,l_dp_qc, l_ad = s01_line2info(tokens)

            ##剩余情况：alt：<*>, 逗号分割的多种情况

            pos, ref_base, alt_base, dp,i16,qs,pl,dp_qc, ad = s02_info_conversion(l_pos, l_ref, l_alt, l_dp,l_i16,l_qs,l_pl,l_dp_qc, l_ad, df_ref)

            if '_' in pos:#仅考虑snp
                continue
            if pos in df_ref_xls.index:
                ###df_ref_xls是否已存在此位置，需判断合并
                new_values = s14_heter(ref_base, alt_base, str(dp),i16,qs,pl,str(dp_qc), ad, df_ref_xls.loc[pos].values,pos)
                df_ref_xls.loc[pos, :] = new_values
            else:
                df_ref_xls.loc[pos, :] = [ref_base, alt_base, str(dp),i16,qs,pl,str(dp_qc), ad]


    vcf.close()
    df_ref_xls['#CHROM']=chrom
    df_ref_xls['POS']=df_ref_xls.index
    df_ref_xls['ID']='.'
    df_ref_xls['QUAL']=0
    df_ref_xls['FILTER']='.'
    df_ref_xls['INFO']=df_ref_xls.apply(lambda x:';'.join(['='.join(['DP',x['DP']]), '='.join(['I16',x['I16']]), '='.join(['QS',x['QS']])]),axis=1)   #DP, I16, QS auxiliary tag for calling
    df_ref_xls['FORMAT']='PL:DP:AD'
    df_ref_xls[sample]=df_ref_xls.apply(lambda x:':'.join([x['PL'],x['DP_qc'],x['AD']]),axis=1)
    df_ref_vcf=df_ref_xls[['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',sample]]
    df_ref_vcf.to_csv(ref_vcf,sep='\t',index=None,header=None)
    return

def run(mpile_vcf,gene_pos_matrix,ref_vcf,chrom,sample):
    mpileupConvertToREF(gene_pos_matrix, mpile_vcf, ref_vcf,chrom,sample)
