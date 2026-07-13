# File-manager thumbnails for .sb models

Stibium's headless renderer can serve as a freedesktop thumbnailer,
so file managers (Nautilus, Thunar, Dolphin via kde-thumbnailer
bridges, ...) show shaded previews of models.

Manual install (packaging integration is a release-milestone task):

    sudo cp stibium-mime.xml /usr/share/mime/packages/
    sudo update-mime-database /usr/share/mime
    sudo cp stibium.thumbnailer /usr/share/thumbnailers/

`stibium` must be on PATH (or edit Exec= to an absolute path).
Clear stale thumbnails with `rm -rf ~/.cache/thumbnails` if previews
don't appear immediately.
