function savefigs(fig, fig_folder, name)

saveas(fig, fig_folder + "\" + name, "epsc");
saveas(gcf, fig_folder + "\" + name, "png");

end