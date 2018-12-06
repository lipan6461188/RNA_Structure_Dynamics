
import tools

tools.scatter_plot(mouse_TS, mouse_TE, color="black", xlabel="TS", ylabel="TE", xlim=[0, 10], ylim=[-2, 2], Parser=mm10_parseTrans, kind='kde')
tools.plt.savefig("figs/mouse_TS_TS_kde.pdf")
tools.plt.show()

