"""Auxiliary review sheet only; official figures remain individual 600 dpi exports."""
import matplotlib.pyplot as plt
def make_overview(root,manifest):
    fig,axs=plt.subplots(4,2,figsize=(12,16));fig.subplots_adjust(wspace=.02,hspace=.08)
    for ax,item in zip(axs.flat,manifest):
        ax.imshow(plt.imread(root/'figures_tables'/item['filename'].split(';')[0]));ax.axis('off')
    axs.flat[-1].axis('off');fig.savefig(root/'qa'/'figure_overview.png',dpi=150);plt.close(fig)
