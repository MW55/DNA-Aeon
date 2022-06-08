import seaborn as sns
import matplotlib.pyplot as plt


def plot_error_prob_for_all(smaller_x: int = 40):
    m_props = {"marker": "x", "markerfacecolor": "black", "markeredgecolor": "black", "markersize": "5"}
    sns.set(style="ticks")
    with open(".OUTFILES/parallel_LT/71766650690_660302228759.fasta", "r") as lt_in:
        lt = lt_in.readlines()[::2]
    with open(".OUTFILES/parallel_Online/63863376951_141614131555.fasta", "r") as online_in:
        online = online_in.readlines()[::2]
    with open(".OUTFILES/parallel_RU10/64174900275_77353324612.fasta", "r") as ru_in:
        ru = ru_in.readlines()[::2]
    tmp = []
    for line in lt:
        val = int(line.split("_")[0].replace(">", ""))
        if val >= smaller_x:
            break
        tmp.append(val)
    f, (ax_box, ax_hist) = plt.subplots(2, sharex=True,
                                        gridspec_kw={"height_ratios": (.15, .85)})
    kwargs = {'cumulative': True, "color": "red"}
    sns.distplot(tmp, hist=True, bins=smaller_x, ax=ax_hist).grid(True)
    ax_hist.set(xlabel="Sequence error probability in %", xlim=(0, smaller_x), ylim=(0, 0.4))
    sns.boxplot(tmp, showmeans=True, meanprops=m_props, ax=ax_box).set_title(
        f"LT (%i Sequences with error probability < %i)" % (len(tmp), smaller_x))
    ax_box.set(yticks=[])
    ax2 = plt.twinx()
    ax2 = sns.distplot(tmp, kde=True, kde_kws=kwargs, hist=False, ax=ax2)
    ax2.set_ylabel('cumulative kde')
    ax2.set(ylim=(0, 1.0))
    sns.despine(ax=ax_hist)
    sns.despine(ax=ax_box, left=True)
    f.show()
    f.savefig(f"05_color_lt_created_error_dist_%i.pdf" % smaller_x, bbox_inches="tight")
    tmp = []
    for line in online:
        val = int(line.split("_")[0].replace(">", ""))
        if val >= smaller_x:
            break
        tmp.append(val)
    f, (ax_box, ax_hist) = plt.subplots(2, sharex=True, gridspec_kw={"height_ratios": (.15, .85)})
    sns.distplot(tmp, hist=True, bins=smaller_x, ax=ax_hist).grid(True)
    ax_hist.set(xlabel="Sequence error probability in %", xlim=(0, smaller_x), ylim=(0, 0.4))
    sns.boxplot(tmp, showmeans=True, meanprops=m_props, ax=ax_box).set_title(
        f"Online (%i Sequences with error probability < %i)" % (len(tmp), smaller_x))
    ax_box.set(yticks=[])
    ax2 = plt.twinx()
    ax2 = sns.distplot(tmp, kde=True, kde_kws=kwargs, hist=False, ax=ax2)
    ax2.set_ylabel('cumulative kde')
    ax2.set(ylim=(0, 1.0))
    sns.despine(ax=ax_hist)
    sns.despine(ax=ax_box, left=True)
    f.show()
    f.savefig(f"05_color_online_created_error_dist_%i.pdf" % smaller_x, bbox_inches="tight")
    tmp = []
    for line in ru:
        val = int(line.split("_")[0].replace(">", ""))
        if val >= smaller_x:
            break
        tmp.append(val)
    f, (ax_box, ax_hist) = plt.subplots(2, sharex=True, gridspec_kw={"height_ratios": (.15, .85)})
    sns.distplot(tmp, hist=True, bins=smaller_x, ax=ax_hist).grid(True)
    ax_hist.set(xlabel="Sequence error probability in %", xlim=(0, smaller_x), ylim=(0, 0.4))
    sns.boxplot(tmp, showmeans=True, meanprops=m_props, ax=ax_box).set_title(
        f"RU10 (%i Sequences with error probability < %i)" % (len(tmp), smaller_x))
    ax_box.set(yticks=[])
    ax2 = plt.twinx()
    ax2 = sns.distplot(tmp, kde=True, kde_kws=kwargs, hist=False, ax=ax2)
    ax2.set_ylabel('cumulative kde')
    ax2.set(ylim=(0, 1.0))
    sns.despine(ax=ax_hist)
    sns.despine(ax=ax_box, left=True)
    f.show()
    f.savefig(f"05_color_ru_created_error_dist_%i.pdf" % smaller_x, bbox_inches="tight")


if __name__ == "__main__":
    plot_error_prob_for_all(100)
