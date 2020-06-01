from ont_fast5_api.fast5_interface import get_fast5_file
from ont_fast5_api.fast5_read import Fast5Read
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import sys

class Fast5ReadWrapper:
    
    """
    1. Fast5ReadWrapper(fast5_filepath, read_id=None, is_load_data=True, 
            basecall_group=None, read_all=False)
        fast5 file need to be basecalled by guppy with flipflop mode. Default
        read the lastest basecalling results (Usually Basecall_1D_000 if you 
        only preform basecalling once). Also can set basecall_group.
        
        If want to read one read, you can set read_id, otherwise it will 
        try to read the first read in fast5_filepath. 
    
        If you want to iter read all read in fast5_filepath, pleas set 
        read_all=True. In this mode, the object will not load data, 
        And then you can use iter_all_read method to iter read each read, 
        each loop return a Fast5ReadWrapper object with a current read_id. 
        And you also can use iter_write_polya to write all polya
        result to a file.
    
    2. load_data()
        extract read data. Usually don't need to use it.
    3. get_read()
        Extract read. Usually don't need to use it. But can use it
        to get_read(read_id) to change read_id (need in fast5_filepath).
        Warning: Need to run close() after processing read.
    4. get_read_ids()
        Extract read_id list in fast5_filepath
    """
    
    CONVERT_BASE_RATIO = 20 # if 0 then not convert
    MIN_POLYA_LENGTH = 15
    MATCH_SCORE = 1
    MISMATCH_SCORE = -1.5
    
    
    def __init__(self,
                 fast5_filepath=None,
                 read_id=None, 
                 read_obj=None,
                 is_load_data=True,
                 basecall_group=None,
                 read_all=False,
                 search_start_base=None, 
                 search_end_base=None,
                 find_base=None,
                 file_search_region=None):        
        self.file = fast5_filepath
        self.basecall_group = None
        self._is_open = False
        self.read_id = read_id
        
        self.search_start_base = search_start_base
        self.search_end_base = search_end_base
        self.find_base = find_base
        self.file_search_region = file_search_region
        self.search_region_df = None
        if search_start_base is None and search_end_base is None:
            if file_search_region:
                self.search_region_df = self.read_search_region(file_search_region)
        
        if not read_all:
            self.load_data()    

    def iter_all_read(self):
        try:
            read_ids = self.read_ids
        except:
            read_ids = self.get_read_ids()
        try:
            for read_id in self.read_ids:
                self.read_id = read_id
                self.load_data(need_close_file=False)
                yield(self)
        finally:
            self.close()
    
    def get_search_region(self):
        if self.search_region_df is not None:
            find_base, search_start_base, search_end_base  = self.search_region_df.loc[self.read_id, ["base", "search_start_base", "search_end_base"]]
            return((find_base, search_start_base, search_end_base))
        else:
            return(self.find_base, self.search_start_base, self.search_end_base)
    
    def read_search_region(self, file_search_region):
        df = pd.read_table(file_search_region)
        df.index = df["read_id"]
        self.search_region_df = df
    
    def iter_write_polya(self, fileout, only_best=True, file_search_region=None):
        if file_search_region:
            self.read_search_region(file_search_region)
            
        dfs = []
        for d in self.iter_all_read():
            df = d.polya_df
            if len(df):
                df["read_id"] = self.read_id
                df["read_length"] = len(self.seq)
                if only_best:
                    df = df.iloc[[df["polya_length"].to_numpy().argmax()]]
                dfs.append(df)
        dfs = pd.concat(dfs)
        dfs = dfs[["read_id", "read_length", "polya_start", "polya_end", 
                            "polya_start_base", "polya_end_base", "polya_length", 
                            "polya_score", "polya_type"]]
        dfs.to_csv(fileout, sep="\t", index=False)
        return(dfs)
            
    def load_data(self, need_get_read=None, need_close_file=True):            
        
        """
        read.get_analysis_attributes(basecall_group) is a dict:
        {'component': 'basecall_1d',
         'model_type': 'flipflop',
         'name': 'ONT Guppy basecalling software.',
         'segmentation': 'Segmentation_000',
         'time_stamp': '2020-03-10T09:56:33Z',
         'version': '3.3.0+ef22818'}
         
        read.get_summary_data(segmentation_name) 
        {'segmentation': {'adapter_max': 0.0,
         'duration_template': 6465,
         'first_sample_template': 766,
         'has_complement': 0,
         'has_template': 1,
         'med_abs_dev_template': 8.23161506652832,
         'median_template': 86.69466400146484,
         'num_events_template': 3232,
         'pt_detect_success': 0,
         'pt_median': 0.0,
        'pt_sd': 0.0}}

        event_table
        ------------------------------------
        |start|base_index|base|event_length|
        |766  |     1    |  G | 12         |
        |778  |     2    |  C | 2          |
        ------------------------------------
        """
        try:
            read = self.get_read(self.read_id)
            if self.basecall_group:
                basecall_group = self.basecall_group
            else:
                basecall_group = read.get_latest_analysis("Basecall_1D")
                
            template_name = basecall_group + "/BaseCalled_template" 
            self.basecall_group = basecall_group
            self.basecall_group_attributes = read.get_analysis_attributes(basecall_group)
            segmentation_name = self.basecall_group_attributes['segmentation']
            if self.basecall_group_attributes['model_type'] != 'flipflop':
                raise ValueError('model type is not flipflop')
            self.raw_data = read.get_raw_data()
            #raw_data: array([805, 496, 514, ..., 521, 531, 643], dtype=int16)
            self.basecall_summary = read.get_summary_data(basecall_group)['basecall_1d_template']
            self.stride = stride = self.basecall_summary['block_stride']
            self.fastq = read.get_analysis_dataset(group_name=template_name, dataset_name='Fastq')
            self.seq = self.fastq.split("\n")[1]
            self.segmentation_summary = read.get_summary_data(segmentation_name)['segmentation']
            self.start = start = self.segmentation_summary['first_sample_template']
            self.move = move = read.get_analysis_dataset(group_name=template_name, dataset_name='Move')
            #move: array([1, 0, 0, ..., 0, 0, 0], dtype=uint8)
            self.event_length = event_length = len(self.move)
            self.end = end = start + (event_length - 1) * stride            
            
            #generate event table
            ls_move_raw_start = (np.arange(event_length)[move==1])*stride + start
            ls_move_raw_start_with_end = np.append(ls_move_raw_start, end)
            ls_event_length = ls_move_raw_start_with_end[1:] - ls_move_raw_start
            self.event_table = pd.DataFrame({"start": ls_move_raw_start, 
                                        "base": list(self.seq), 
                                        "event_length": ls_event_length})
            self.samples_per_nt = np.mean(ls_event_length[ls_event_length <=np.quantile(ls_event_length, 0.95)])
            self.find_polyA_or_polyT()
            #to numpy to increase speed (for raw2baseindex)
            #self.ls_move_raw_start_with_end = ls_move_raw_start_with_end
            #self.ls_move_raw_start = ls_move_raw_start
            #self.base_index = np.arange(ls_move_raw_start)
        finally:
            if need_close_file:
                self.close()
    
    #def raw2baseindex(self, start):
    #    return(self.base_index[self.ls_move_raw_start == start])
    
    def savefig(self, filename):
        self.fig.savefig(filename)
    
    def plot(self, figsize = None, plot_all_polyA_region=False, plot_base=False, plot_base_line=False, xlim=None, ylim=None):
        
        raw_data = self.raw_data
        event_table = self.event_table #don't change
        end = self.end
        start = self.start
        
        plot_polya = True
        try:
            polya_start = self.polya_start
            polya_end = self.polya_end
            polya_type = self.polya_type
            polya_length = self.polya_length
            polya_df = self.polya_df
        except:
            plot_polya = False

        if figsize:
            fig, ax = plt.subplots(figsize=figszie)
        else:
            fig, ax = plt.subplots()
            
        ax.plot(np.arange(len(raw_data)), raw_data)
        ax.axvspan(start, end, facecolor='g', alpha=0.25)
        if plot_polya:
            ax.set_title("%s %s-%s length: %.1f" % (polya_type, polya_start, polya_end, polya_length))
            if plot_all_polyA_region:
                for i, d in polya_df.iterrows():
                    if d["polya_type"] == "polyA":
                        color = "r"
                    else:
                        color = "b"
                    ax.axvspan(d["polya_start"], d["polya_end"], facecolor=color, alpha=0.25)
            else:
                ax.axvspan(polya_start, polya_end, facecolor='r', alpha=0.25)                    

        if plot_base:
            starts = event_table["start"]
            ends = np.append(event_table["start"][1:], end)
            bases = event_table["base"]
            for start, end, base in zip(starts, ends, bases):
                y_pos = raw_data[start:end].max() + 10
                ax.text((start+end-1)/2, y_pos, base, clip_on=True, horizontalalignment="center",
                       verticalalignment='center') #base + str(end-start)
                if plot_base_line:
                    ax.axvline(start, color="grey", zorder=1)
            if plot_base_line:
                    ax.axvline(end, color="grey", zorder=1)
        
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylims(ylim)
            
        self.ax = ax
        self.fig = fig
    
    def find_polyA_or_polyT(self):
        
        find_base, start_base, end_base = self.get_search_region()
        
        e = self.event_table
        samples_per_nt = self.samples_per_nt
        
        if find_base:
            result, df = self.find_polyA_by_event_data(find_base, start_base, end_base)
        else:            
            polya_result, polya_df = self.find_polyA_by_event_data("A", start_base, end_base)
            polyt_result, polyt_df = self.find_polyA_by_event_data("T", start_base, end_base)
            result = polya_result if polya_result[-3] >= polyt_result[-3] else polyt_result
            df = pd.concat([polya_df, polyt_df])
        (self.polya_start, self.polya_end, self.polya_start_base, 
             self.polya_end_base, self.polya_length, self.polya_score,
             self.polya_type) = result
        
        self.polya_df = df
        return (result, df)
    
    def find_polyA_by_event_data(self, find_base="A", 
                                 start_base=None, 
                                 end_base=None):
    
        """
        e is not changed.
        
        convert base with event_length >= convert_base_ratio*samples_per_nt and 
        adjecent with find_base to find_base. if convert_base_ratio == 0, then not convert.
        use numpy but not series to increase speed.
        
        results:
        [result, df]
        result is the max polya length region (even the polya length < MIN_POLYA_LENGTH):
        [polya_start, polya_end, polya_start_base, polya_end_base, polya_length, polya_score, polya_type]
        polya_start, polya_end: raw signal start of base, 0-based
        polya_start_base, polya_end_base: 1-based base index (not base A T C G, but position)
        polyA_type: "polyA" or "polyT"
        df is DataFrame storing polya regions with polya length more than MIN_POLYA_LENGTH.
        df columns:
        "polya_start", "polya_end", "polya_start_base", 
        "polya_end_base", "polya_length", "polya_score", "polya_type"
        """        
        
        e = self.event_table
        samples_per_nt = self.samples_per_nt
        convert_base_ratio = self.CONVERT_BASE_RATIO
        min_polya_length = self.MIN_POLYA_LENGTH
        match_score = self.MATCH_SCORE
        mismatch_score = self.MISMATCH_SCORE
        
        df = pd.DataFrame(columns=["polya_start", "polya_end", 
                            "polya_start_base", "polya_end_base", "polya_length", 
                            "polya_score", "polya_type"])
        
        base_index = np.arange(len(e))
        
        #if fix search region, select and copy data.
        if start_base:
            if not end_base: end_base = len(e)
            e = e[(start_base-1):end_base].copy()
            base_index = base_index[(start_base-1):end_base] # auto copy
        else:
            if end_base:
                e = e[:end_base].copy()
                base_index = base_index[:end_base]
        last_base_index = base_index[-1]
        
        end = e.iloc[-1, :]["start"] + e.iloc[-1, :]["event_length"]
        base_is_right = e["base"].to_numpy() == find_base
        
        #convert adjecent long base
        if convert_base_ratio:
            f = np.logical_and(e["event_length"] >= samples_per_nt*convert_base_ratio,
                   np.logical_not(base_is_right))
            if f.sum() >0:
                for i in np.where(f):
                    try:
                        if base_is_right[i-1]: 
                            base_is_right[i] = True
                            next
                    except:
                        pass
                    try:
                        if base_is_right[i+1]:
                            base_is_right[i] = True
                    except:
                        pass
                
        #group adjecent same bases
        f = base_is_right != np.insert(base_is_right[:-1], 0, True)
        f[0] = True

        base_is_right = base_is_right[f]
        total_seed_num = base_is_right.sum()
        if total_seed_num > 0:
            starts = e["start"][f].to_numpy()
            starts_with_end = np.append(starts, end)
            event_lengths = starts_with_end[1:] - starts
            base_scores = event_lengths * np.where(base_is_right, match_score, mismatch_score)
            base_index = base_index[f]

            seed_scores = np.zeros((total_seed_num, 3), dtype=int)

            if base_is_right[0]:
                init_i = 0
            else:
                init_i = 1

            now_s, now_score, now_seed_num = init_i, base_scores[init_i], 0
            for i in np.arange(init_i, len(base_scores)-2,2):
                #not greedy search
                #i is match(A), i+1 is mismatch, i+2 is match
                if (now_score + base_scores[i+1] > 0) and (base_scores[i+2] + base_scores[i+1] > 0):
                    now_score = now_score + base_scores[i+2] + base_scores[i+1]
                else:
                    seed_scores[now_seed_num] = (now_s, i, now_score)
                    now_s = i + 2
                    now_score = base_scores[i+2]
                    now_seed_num += 1
            seed_scores[now_seed_num] = (now_s, i+2, now_score)

            seed_event_starts = starts_with_end[seed_scores[:(now_seed_num+1),0]]
            seed_event_ends = starts_with_end[seed_scores[:(now_seed_num+1),1] + 1]
            seed_event_length = seed_event_ends - seed_event_starts
            seed_polya_length = seed_event_length/samples_per_nt
                        
            f = np.where(seed_polya_length > min_polya_length)
            df = pd.DataFrame(seed_scores[f], columns=["s", "e", "polya_score"])
            df["polya_start"] = seed_event_starts[f]
            df["polya_end"] = seed_event_ends[f] - 1 # need minus 1  debug
            df["polya_start_base"] = base_index[df["s"]]+1
            df["polya_end_base"] = np.append(base_index, last_base_index)[df["e"] + 1]
            df["polya_length"] = seed_polya_length[f]
            df["polya_type"] = "poly" + find_base
            df = df[["polya_start", "polya_end", "polya_start_base", "polya_end_base", "polya_length", "polya_score", "polya_type"]]
                        
            i = seed_event_length.argmax()
            start, end, score = seed_scores[i, :]
            result = (seed_event_starts[i], seed_event_ends[i], 
                    base_index[start]+1, base_index[end]+1,
                    seed_event_length[i]/samples_per_nt, score, "poly"+find_base)
            
        else:
            result = (0, 0, 0, 0, 0, 0, "poly"+find_base)
        
        if not len(df.index):
            df = pd.DataFrame([result], columns=["polya_start", "polya_end", 
                        "polya_start_base", "polya_end_base", "polya_length", 
                        "polya_score", "polya_type"])
        #if return_only_best:
        #        df = df.iloc[[df["polya_length"].to_numpy().argmax()]]
                        
        self.polya_df = df
        self.polya_result = result
        
        (self.polya_start, self.polya_end, self.polya_start_base, 
             self.polya_end_base, self.polya_length, self.polya_score,
             self.polya_type) = result
        
        return(result, df)

    def get_read_ids(self, need_to_close=True):
        self.open()
        try:
            read_ids = self.io.get_read_ids()
            self.read_ids = read_ids
        finally:
            if need_to_close:
                self.close()
        return(read_ids)
    
    def get_read(self, read_id=None):
        self.open()
        if not read_id:
            read_id = self.read_id
        if not read_id:
            read_ids = self.get_read_ids(need_to_close=False)
            try:
                read_id = read_ids[0]
            except:
                pass
        self.read_id = read_id
        if read_id:
            read = self.io.get_read(read_id)
        else:
            read = None
        return(read)
            
    def open(self, force_open=False):
        if force_open or (not self._is_open):
            try:
                self.close()
            except:
                pass
            self.io = get_fast5_file(self.file, mode="r")
        self._is_open = True
        
    def close(self):
        self._is_open = False
        try:
            self.io.close()
        except:
            pass
        
    def find_polyA(self):
        pass
        
def extract_fast5_polyA(filein, fileout, only_best=True, file_search_region=None):
    Fast5ReadWrapper(filein, read_all=True).iter_write_polya(fileout, only_best, file_search_region)
    
    
def main(argv):
    filein, fileout = sys.argv[1:3]
    try:
        only_best = int(sys.argv[3])
    except:
        only_best = 1
    try:
        file_search_region = sys.argv[4]
    except:
        file_search_region = None
    extract_fast5_polyA(filein, fileout, only_best, file_search_region)
    
if __name__ == "__main__":
    main(sys.argv)
