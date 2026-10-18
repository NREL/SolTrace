from datetime import datetime
import time, pickle
from dataclasses import dataclass
import numpy as np
import orjson
from colorama import just_fix_windows_console, Fore, Back, Style
just_fix_windows_console()

right_float = lambda fl, buf, p: f"{f'{fl:.{p}f}':>{buf}}"
right_int = lambda i, buf: f"{str(i):>{buf}}"
right_str = lambda s, buf: f"{s:>{buf}}"
left_float = lambda fl, buf, p: f"{f'{fl:.{p}f}':<{buf}}"

def _min_max_n_values(arr, n):
    sort = np.sort(arr)
    return sort[:n], sort[-n:]

def _fmt_header(title: str, cols: list[str], key_buf: int, col_buf: int):
    lines = [f'{title:<{key_buf}}    ' + ''.join(f'{c:>{col_buf}}' for c in cols)]
    lines.append('-'*(len(lines[0])))
    return lines

def _highlight_float(val: float, buf: int, p: int, color_bounds: tuple[list[float], list[float], float, float]):
    if val == color_bounds[2]:
        return f'{Fore.LIGHTGREEN_EX}{right_float(val, buf, p)}{Style.RESET_ALL}'
    elif val == color_bounds[3]:
        return f'{Fore.RED}{right_float(val, buf, p)}{Style.RESET_ALL}'
    elif val in color_bounds[0]:
        return f'{Fore.YELLOW}{right_float(val, buf, p)}{Style.RESET_ALL}'
    elif val in color_bounds[1]:
        return f'{Fore.MAGENTA}{right_float(val, buf, p)}{Style.RESET_ALL}'
    else: return right_float(val, buf, p)

def _highlight_int(val: int, buf: int, color_bounds: tuple[list[float], list[float], float, float]):
    if val == color_bounds[2]:
        return f'{Fore.LIGHTGREEN_EX}{right_int(val, buf)}{Style.RESET_ALL}'
    elif val == color_bounds[3]:
        return f'{Fore.RED}{right_int(val, buf)}{Style.RESET_ALL}'
    elif val in color_bounds[0]:
        return f'{Fore.YELLOW}{right_int(val, buf)}{Style.RESET_ALL}'
    elif val in color_bounds[1]:
        return f'{Fore.MAGENTA}{right_int(val, buf)}{Style.RESET_ALL}'
    else: return right_int(val, buf)

def _highlight(toggle: bool, val: float | int, buf: int, p: int = 7, color_bounds: list[tuple[list[float], list[float], float, float]] = None, color_i: int = 0):
    is_float = isinstance(val, float)
    if toggle:
        buf_func = _highlight_float if is_float else _highlight_int
        _args = (val, buf, p, color_bounds[color_i]) if is_float else (val, buf, color_bounds[color_i])
    else:
        buf_func = right_float if is_float else right_int
        _args = (val, buf, p) if is_float else (val, buf)
    return buf_func(*_args)

class timer():
    def __init__(self, verbose: bool = False, precision: int = 7):
        self.start = {}
        self.verbose = verbose
        self.history = { 'total': 0 }
        self.precision = precision
    
    # amdahl util

    def __repr__(self):
        if not len(self.history) > 1: return 'Timer has no times to report.'

        # get info to print
        timed_keys, stats, stats_cols = self.summary()
        # toggle at over 12 so fastest/slowest 3 are highlighted
        color_toggle = stats.shape[0] > 11
        colored_values = []

        if color_toggle:
            threshold = int(np.ceil(.2 * stats.shape[0]))
            for i in range(stats.shape[1]):
                _mins, _maxs = _min_max_n_values(stats[:, i], threshold)
                colored_values.append((_mins, _maxs, _mins[0], _maxs[-1]))

        # find the buffers for keys and total col to align the values properly
        max_key_len = max(len(str(k)) for k in self.history.keys())
        max_sum_len = max(len(f'{_sum:.{self.precision}f}') for _sum in stats[:, 0]) + 2

        # format header
        formatted_lines = _fmt_header('timer summary:', stats_cols, max_key_len, max_sum_len)

        # for each row in the stats array add key name and timing stats
        # if enough events timed highlight values
        for i in range(stats.shape[0]):
            vals = stats[i, :]
            report_line = [f"  {timed_keys[i]:<{max_key_len}}  "]
            report_line.append(_highlight(color_toggle, vals[0], max_sum_len, self.precision, colored_values, 0)) # [0] if color_toggle else None
            if (count := int(vals[-1])) > 1:
                report_line.extend([_highlight(color_toggle, vals[j], max_sum_len, self.precision, colored_values, j)
                                    for j in range(1, stats.shape[1] - 1)])
                report_line.append(_highlight(color_toggle, count, max_sum_len, self.precision, colored_values, -1))
            formatted_lines.append(''.join(report_line))

        # add total row and footer
        report_str = f'{right_float(self.history['total'] / 60, max_sum_len)} [m]' if self.history['total'] > 60 else f'{right_float(self.history['total'], max_sum_len, self.precision)} [s]'
        formatted_lines.append(f"  {'total':<{max_key_len}}  {report_str}")
        formatted_lines.append(formatted_lines[1])
        
        return "\n".join(formatted_lines)
    
    def __iadd__(self, other: dict):
        # if not a dictionary skip
        if type(other) != dict: return self
        # doesn't add to total because it's assumed that this info is coming 
        # from somewhere timer can't be called, ie in an executable, and it
        # is assumed that the executable call is wrapped in t.ic t.oc
        for key, diff in other.items():
            s = f'  - {key}'
            if s in self.history: self.history[s].append(diff)
            else:                 self.history[s] = [diff]
        return self

    def ic(self, s: str = ''):
        self.start[s] = time.perf_counter()
        if self.verbose: print(f'started{" " + s if s else ""} at: {datetime.now().strftime("%Y/%m/%d %H:%M:%S")}')

    def oc(self, s: str = ''):
        assert s in self.start, f"{s} never started"

        diff = time.perf_counter() - self.start[s]
        if self.verbose: print(f'ended{" " + s if s else ""} at: {datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n{s + " " if s else ""}took: {diff}')
        if s: 
            self.history['total'] += diff
            if s in self.history: self.history[s].append(diff)
            else:                 self.history[s] = [diff]
        return diff

    def summarize(self, s: str):
        assert s in self.history, f"{s} never timed"

        diffs = self.history[s]
        value = sum(diffs)
        c     = len(diffs)

        ave    = None
        stddev = None
        median = None
        if (c > 1):
            ave    = value / c
            stddev = np.std(diffs)
            median = np.median(diffs)

        return value, ave, stddev, median, c

    def summary(self):
        timed_keys = [k for k in self.history.keys() if k != 'total']
        stats = np.array([self.summarize(k) for k in timed_keys])
        stats_cols = ["total", "average", "std", "median", "counts"]
        return timed_keys, stats, stats_cols
    
    def to_json(self, fname: str):
        f = open(fname, mode='wb')
        f.write(orjson.dumps(self.history, option=orjson.OPT_INDENT_2))
        f.close()

    def load_json(self, fname: str):
        f = open(fname, mode='rb')
        self.history = orjson.loads(f.read())
        f.close()

@dataclass
class benchmark_record:
    max:    float | int = -np.inf # historical maximum
    min:    float | int = np.inf  # historical minimum
    recent: float | int = np.nan  # most recently recorded
    time:   datetime    = None    # time of recent was recorded

    def __repr__(self):
        return str(self.recent)

    def __format__(self, format_spec):
        if not format_spec: format_spec = '.2f'
        if isinstance(self.recent, int): format_spec = ''
        fmt = '{:' + format_spec + '}'
        return fmt.format(self.recent)

    def update(self, val: float | int, t: datetime):
        self.recent = val
        self.time = t
        if val > self.max: self.max = val
        if val < self.min: self.min = val

    def highlight(self, val: float | int, buf: int, p: int = 7):
        buf_func = right_float if isinstance(self.recent, float) else right_int
        _args = (val, buf, p) if isinstance(self.recent, float) else (val, buf)

        if val < self.min:
            return f'{Fore.LIGHTGREEN_EX}{buf_func(*_args)}{Style.RESET_ALL}'
        elif val > self.max:
            return f'{Fore.RED}{buf_func(*_args)}{Style.RESET_ALL}'
        elif val < self.recent:
            return f'{Fore.YELLOW}{buf_func(*_args)}{Style.RESET_ALL}'
        elif val > self.recent:
            return f'{Fore.MAGENTA}{buf_func(*_args)}{Style.RESET_ALL}'
        else: return buf_func(*_args)

class benchmark_store:
    def __init__(self):
        self.store: dict[str, dict[str, benchmark_record]] = {}
        self.benchmarks: list[str] = []
        self.stats_cols: list[str] = []

    def create(self, benchmarks: list[str], stats: np.ndarray, stats_cols: list[str], t: datetime):
        assert len(benchmarks) == stats.shape[0], f"benckmarks and stats provided are different lengths ({len(benchmarks)} != {stats.shape[0]})"
        assert len(stats_cols) == stats.shape[1], f"stats and stats names provided are different lengths ({len(stats_cols)} != {stats.shape[1]})"

        self.store: dict[str, dict[str, benchmark_record]] = {}
        self.benchmarks: list[str] = benchmarks
        self.stats_cols: list[str] = stats_cols

        for i in range(stats.shape[0]):
            self.store[benchmarks[i]] = {}
            for j in range(stats.shape[1]):
                self.store[benchmarks[i]][stats_cols[j]] = benchmark_record(stats[i, j], stats[i, j], stats[i, j], t)
        return self

    def dump(self, filename: str):
        f = open(filename, mode='wb')
        _dump = {
            'store':      self.store,
            'benchmarks': self.benchmarks,
            'stats_cols': self.stats_cols,
        }
        pickle.dump(_dump, f)
        f.close()

    def load(self, filename: str):
        f = open(filename, mode='rb')
        _load = pickle.load(f)
        f.close()

        assert isinstance(_load, dict), 'expects a dictionary'
        assert 'store' in _load,        'expect store to be in loaded dictionary'
        assert 'benchmarks' in _load,   'expect benchmarks to be in loaded dictionary'
        assert 'stats_cols' in _load,   'expect stats_cols to be in loaded dictionary'

        self.store      = _load['store']
        self.benchmarks = _load['benchmarks']
        self.stats_cols = _load['stats_cols']

        return self

    def update(self, benchmarks: list[str], stats: np.ndarray, stats_cols: list[str], t: datetime):
        assert len(benchmarks) == stats.shape[0], f"benckmarks and stats provided are different lengths ({len(benchmarks)} != {stats.shape[0]})"
        assert len(stats_cols) == stats.shape[1], f"stats and stats names provided are different lengths ({len(stats_cols)} != {stats.shape[1]})"

        self.__check_append_keys(stats_cols)        

        for i in range(stats.shape[0]):
            benchmark = benchmarks[i]

            # init inner dictionary in necessary
            if benchmark not in self.benchmarks:
                self.store[benchmark] = { c: benchmark_record() for c in self.stats_cols }
                self.benchmarks.append(benchmark)

            # update benchmarks stats
            for j in range(stats.shape[1]):
                self.store[benchmark][self.stats_cols[j]].update(stats[i, j], t)

    def __check_append_keys(self, stats_cols: list[str]):
        # find any possible new stats to store
        possible_new_cols = []
        for col in stats_cols:
            if col not in self.stats_cols:
                possible_new_cols.append(col)

        # if no new cols don't do anything else
        if len(possible_new_cols) == 0: return

        # add new cols to list and add new records to store
        self.stats_cols.extend(possible_new_cols)
        for new_col in possible_new_cols:
            for k in self.store.keys():
                self.store[k][new_col] = benchmark_record()

    def compare(self, benchmarks: list[str], stats: np.ndarray, stats_cols: list[str], p: int = 7) -> str:
        common_benchmarks = sorted(list(set(benchmarks) & set(self.benchmarks)))
        common_stats_cols = sorted(list(set(stats_cols) & set(self.stats_cols)))

        _is = { _b: benchmarks.index(_b) for _b in common_benchmarks }
        _js = { _c: stats_cols.index(_c) for _c in common_stats_cols }

        # find the buffers for keys and total col to align the values properly
        max_key_len = max(len(str(_b)) for _b in common_benchmarks)
        max_sum_len = max(len(f'{_sum:.{p}f}') for _sum in stats[:, _js[common_stats_cols[0]]]) + 2

        # format header
        formatted_lines = _fmt_header('store comparison:', common_stats_cols, max_key_len, 2 * max_sum_len + 1)

        for _b in common_benchmarks:
            report_line = [f"  {_b:<{max_key_len}}  "]

            for _c in common_stats_cols:
                _rec = self.store[_b][_c]
                report_line.append(f'{_rec.highlight(stats[_is[_b], _js[_c]], max_sum_len, p)} ({_rec:.{p}f})')
            formatted_lines.append(''.join(report_line))

        formatted_lines.append(formatted_lines[1])
        return "\n".join(formatted_lines)

