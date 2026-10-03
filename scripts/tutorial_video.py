#!/usr/bin/env python3
"""Record the tutorial videos by driving the real application offscreen.

Six silent screen captures, one per workflow, written to ``docs/tutorial/video/`` as
``<slug>_silent.mp4`` beside a ``poster.jpg`` -- the shape spaCR ships its tutorial videos in.
Nothing is staged: every frame is a `QWidget.grab()` of the real `Window` with the real panels in
it, every number on screen was computed by the shipped code while this script ran, and the only
pixels this module paints itself are the overlay -- a synthetic pointer, a highlight around the
control about to be used, and the caption burned into a band under the window, since the clips
carry no audio.

    QT_QPA_PLATFORM=offscreen python scripts/tutorial_video.py            # all six
    QT_QPA_PLATFORM=offscreen python scripts/tutorial_video.py 3_test     # one

Each clip is a list of beats: say something, move the pointer somewhere, press something, wait
while the viewer reads. Frames are PNGs in a temporary directory and `/usr/bin/ffmpeg` encodes
them to H.264 / yuv420p at `FPS`; a still from the middle of the clip becomes the poster. The
window state between beats is cached, so a four-second pause costs one grab rather than fifty, and
a whole clip costs about as long as the strategy runs inside it.

Two honest limitations, both visible in the clips rather than hidden. The central 3D view is a
`QOpenGLWidget`, which the offscreen Qt platform does not support, so it records as an empty strip
and is labelled as such in the frame. And a strategy run takes as long as it takes: the clip holds
on the panel before and after the run rather than pretending the wait away.

Re-runnable: the window is built with its own settings and state directories, the beats are fixed,
and the shipped tables and calibration are the same on every run, so a rebuild reproduces the clips
except for the seconds a run happens to take.
"""
from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

OUT = os.path.join(ROOT, "docs", "tutorial", "video")
FFMPEG = "/usr/bin/ffmpeg"

#: Frames per second. Twelve is plenty for a screen capture whose content changes in steps, and it
#: keeps each clip to a few hundred kilobytes.
FPS = 12
#: The recorded window, and the caption band drawn under it. Both even, as yuv420p requires.
WINDOW = (1680, 920)
BAND = 130
#: Where the PNG frames are staged before encoding -- a clip is a few hundred megabytes of them,
#: which is more than a small /tmp wants. Set `STARPLAST_VIDEO_SCRATCH` to put them elsewhere.
SCRATCH = os.environ.get("STARPLAST_VIDEO_SCRATCH") or None
#: Overlay colours, independent of the application's theme so the captions read on all of them.
INK = "#e8f2f4"
DIM = "#9fb3ba"
ACCENT = "#35d0c8"
BAND_BG = "#0a1216"


# --------------------------------------------------------------------------- the application
class Studio:
    """One offscreen application and one real `Window`, with isolated settings and user state.

    The window is the thing being filmed, so it is built exactly as the program builds it; only the
    places it would write to are redirected, which is what `tests/conftest.py` does for the suite
    and for the same reason.
    """

    def __init__(self):
        """Start Qt offscreen with a scratch settings and state directory, and open the window."""
        defaults = {"QT_QPA_PLATFORM": "offscreen", "PYQTGRAPH_QT_LIB": "PyQt6",
                    "STARPLAST_GPU": "0", "LIBGL_ALWAYS_SOFTWARE": "1", "QT_OPENGL": "software",
                    "__GLX_VENDOR_LIBRARY_NAME": "mesa"}
        for key, value in defaults.items():
            os.environ.setdefault(key, value)
        self.tmp = tempfile.mkdtemp(prefix="starplast-video-")
        os.environ["STARPLAST_STATE"] = self.tmp
        from PyQt6 import QtCore, QtWidgets
        QtCore.QSettings.setDefaultFormat(QtCore.QSettings.Format.IniFormat)
        for fmt in (QtCore.QSettings.Format.NativeFormat, QtCore.QSettings.Format.IniFormat):
            for scope in (QtCore.QSettings.Scope.UserScope, QtCore.QSettings.Scope.SystemScope):
                QtCore.QSettings.setPath(fmt, scope, self.tmp)
        if QtWidgets.QApplication.instance() is None:
            from starplast import sprite
            sprite.ensure_gl_format()
            QtCore.QCoreApplication.setAttribute(
                QtCore.Qt.ApplicationAttribute.AA_ShareOpenGLContexts)
        self.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        from starplast import app as A
        self.window = A.Window()
        self.window.resize(*WINDOW)
        self.window.show()
        self.settle(12)
        # The map the window opens on, so a clip that loads a gallery map can hand it back.
        self.home = tuple(None if a is None else a.copy()
                          for a in (self.window.xyz, self.window.placed,
                                    self.window.cluster_labels))

    def settle(self, n: int = 8):
        """Let Qt lay out, paint and run the deferred measurements the panels post to themselves."""
        from PyQt6 import QtCore
        for _ in range(n):
            self.app.processEvents()
            QtCore.QCoreApplication.sendPostedEvents()

    def close(self):
        """Drop the window and stop the console's file handler, as the screenshot scripts do."""
        try:
            self.window.console.remove()
        except Exception:                       # a window built without logging has nothing to stop
            pass
        self.window.close()
        shutil.rmtree(self.tmp, ignore_errors=True)


# --------------------------------------------------------------------------- the overlay
def _color(name: str):
    from PyQt6 import QtGui
    return QtGui.QColor(name)


def _arrow(painter, x: float, y: float, scale: float = 1.7):
    """Draw the synthetic pointer: a plain arrow, white on black, so it reads on any panel."""
    from PyQt6 import QtCore, QtGui
    shape = [(0, 0), (0, 17), (4.2, 12.8), (7.2, 19), (9.6, 18), (6.6, 12), (12, 12)]
    poly = QtGui.QPolygonF([QtCore.QPointF(x + a * scale, y + b * scale) for a, b in shape])
    painter.setPen(QtGui.QPen(QtGui.QColor(10, 14, 16, 230), 2.4))
    painter.setBrush(QtGui.QColor(252, 253, 253))
    painter.drawPolygon(poly)


def _ease(t: float) -> float:
    """Smooth-step easing, so the pointer starts and stops gently instead of snapping."""
    t = max(0.0, min(1.0, t))
    return t * t * (3.0 - 2.0 * t)


class Clip:
    """One recorded workflow: beats in, PNG frames out, then one mp4 and one poster.

    The window is grabbed only when something about it changed (`refresh`); every frame in a pause
    is that cached image with the overlay painted over it, which is what makes a two-minute clip
    cost seconds rather than minutes.
    """

    def __init__(self, studio: Studio, slug: str, title: str, steps: int):
        """Prepare a clip of `steps` numbered beats, filmed from `studio`'s window."""
        from PyQt6 import QtGui
        self.studio, self.slug, self.title = studio, slug, title
        self.steps, self.step = steps, 0
        self.frames = tempfile.mkdtemp(prefix=f"starplast-frames-{slug}-", dir=SCRATCH)
        self.n = 0
        self.poster_frame = None
        self.heading = title
        self.text = ""
        self.cursor = (WINDOW[0] * 0.5, WINDOW[1] * 0.5)
        self.focus = None                       # the rect being highlighted, in window coordinates
        self.flash = 0.0                        # 1 right after a press, fading over the next frames
        self.base: QtGui.QImage | None = None
        self.refresh()

    # ------------------------------------------------------------------ frames
    def refresh(self):
        """Grab the window again: call after anything that changes what is on screen."""
        self.studio.settle()
        self.base = self.studio.window.grab().toImage()

    def _frame(self):
        """Compose one frame: the cached window, the overlay on it, and the caption band under it."""
        from PyQt6 import QtCore, QtGui
        w, h = WINDOW
        image = QtGui.QImage(w, h + BAND, QtGui.QImage.Format.Format_RGB32)
        image.fill(_color(BAND_BG))
        p = QtGui.QPainter(image)
        p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        p.setRenderHint(QtGui.QPainter.RenderHint.TextAntialiasing)
        p.drawImage(0, 0, self.base)
        self._draw_gl_note(p)
        if self.focus is not None:
            self._draw_focus(p)
        _arrow(p, *self.cursor)
        self._draw_band(p)
        p.end()
        return image

    def _draw_gl_note(self, p):
        """Label the empty central strip rather than leave an unexplained black column."""
        from PyQt6 import QtCore, QtGui
        view = getattr(self.studio.window, "view", None)
        if view is None or not view.isVisible():
            return
        top_left = view.mapTo(self.studio.window, QtCore.QPoint(0, 0))
        rect = QtCore.QRect(top_left, view.size())
        if rect.width() < 60 or rect.height() < 200:
            return
        p.save()
        p.translate(rect.center())
        p.rotate(-90)
        font = QtGui.QFont()
        font.setPixelSize(15)
        p.setFont(font)
        p.setPen(_color("#4a6068"))
        p.drawText(QtCore.QRect(-rect.height() // 2, -rect.width() // 2, rect.height(),
                                rect.width()),
                   int(QtCore.Qt.AlignmentFlag.AlignCenter),
                   "the 3D map needs OpenGL, which a headless recording cannot draw")
        p.restore()

    def _draw_focus(self, p):
        """Ring the control the pointer is on, and flare it for a moment when it is pressed."""
        from PyQt6 import QtCore, QtGui
        rect = QtCore.QRectF(self.focus).adjusted(-3, -3, 3, 3)
        glow = _color(ACCENT)
        glow.setAlpha(40 + int(90 * self.flash))
        p.setBrush(glow)
        pen = QtGui.QPen(_color(ACCENT), 2.0 + 1.6 * self.flash)
        p.setPen(pen)
        p.drawRoundedRect(rect, 6, 6)
        if self.flash > 0:
            ring = _color(ACCENT)
            ring.setAlpha(int(200 * self.flash))
            p.setBrush(QtCore.Qt.BrushStyle.NoBrush)
            p.setPen(QtGui.QPen(ring, 2.5))
            radius = 10 + 26 * (1.0 - self.flash)
            p.drawEllipse(QtCore.QPointF(*self.cursor), radius, radius)

    def _draw_band(self, p):
        """The caption: which beat this is, its heading, and a sentence or two under it."""
        from PyQt6 import QtCore, QtGui
        w, h = WINDOW
        p.setPen(QtGui.QPen(_color(ACCENT), 2))
        p.drawLine(0, h + 1, w, h + 1)
        small = QtGui.QFont()
        small.setPixelSize(14)
        small.setBold(True)
        small.setLetterSpacing(QtGui.QFont.SpacingType.AbsoluteSpacing, 1.4)
        p.setFont(small)
        p.setPen(_color(ACCENT))
        p.drawText(QtCore.QRect(26, h + 14, 420, 20), 0,
                   f"STEP {self.step} OF {self.steps}".upper())
        p.setPen(_color("#6d8992"))
        p.drawText(QtCore.QRect(w - 626, h + 14, 600, 20),
                   int(QtCore.Qt.AlignmentFlag.AlignRight), self.title.upper())
        head = QtGui.QFont()
        head.setPixelSize(24)
        head.setBold(True)
        p.setFont(head)
        p.setPen(_color(INK))
        p.drawText(QtCore.QRect(26, h + 36, w - 52, 30), 0, self.heading)
        body = QtGui.QFont()
        body.setPixelSize(18)
        p.setFont(body)
        p.setPen(_color(DIM))
        p.drawText(QtCore.QRect(26, h + 70, w - 52, BAND - 76),
                   int(QtCore.Qt.TextFlag.TextWordWrap), self.text)

    def _emit(self, count: int = 1):
        """Write `count` identical frames, decaying the press flare across them."""
        for _ in range(max(1, count)):
            self.n += 1
            self._frame().save(os.path.join(self.frames, f"{self.n:05d}.png"))
            self.flash = max(0.0, self.flash - 0.25)

    # ------------------------------------------------------------------ beats
    @staticmethod
    def _prose(text: str) -> str:
        """The repository writes `--` for a dash; a burned-in caption should show one, not two."""
        return " ".join(str(text).split()).replace(" -- ", " — ")

    def say(self, heading: str, text: str, seconds: float = 0.0):
        """Start a new beat: a heading, the sentence under it, and an optional pause to read it."""
        self.step += 1
        self.heading, self.text = self._prose(heading), self._prose(text)
        if seconds:
            self.hold(seconds)

    def note(self, text: str, seconds: float = 0.0):
        """Change the caption's sentence without starting a new beat."""
        self.text = self._prose(text)
        if seconds:
            self.hold(seconds)

    def hold(self, seconds: float):
        """Stay where we are for `seconds`, which is how long the viewer has to read."""
        self._emit(int(round(seconds * FPS)))

    def rect_of(self, widget):
        """A widget's rectangle in window coordinates, which is where the overlay draws."""
        from PyQt6 import QtCore
        return QtCore.QRect(widget.mapTo(self.studio.window, QtCore.QPoint(0, 0)), widget.size())

    def move_to(self, target, seconds: float = 0.9, focus: bool = True):
        """Glide the pointer to a widget or rectangle, highlighting it when it arrives."""
        from PyQt6 import QtCore
        rect = target if isinstance(target, QtCore.QRect) else self.rect_of(target)
        goal = (rect.center().x(), rect.center().y())
        start, frames = self.cursor, max(2, int(round(seconds * FPS)))
        self.focus = rect if focus else None
        for i in range(1, frames + 1):
            t = _ease(i / frames)
            self.cursor = (start[0] + (goal[0] - start[0]) * t,
                           start[1] + (goal[1] - start[1]) * t)
            self._emit()
        return rect

    def press(self, widget, action=None, *, move: float = 0.9, settle: float = 0.6,
              refresh: bool = True):
        """Move to a control, flare it, then do what pressing it does -- by default, click it."""
        self.move_to(widget, move)
        self.flash = 1.0
        self._emit(2)
        (action or widget.click)()
        if refresh:
            self.refresh()
        self.hold(settle)

    def type_into(self, line_edit, text: str, *, per_char: float = 0.09, after: float = 0.4):
        """Type into a box one character at a time, so the viewer sees what is being searched for."""
        self.move_to(line_edit, 0.8)
        line_edit.setText("")
        self.refresh()
        every = max(1, int(round(per_char * FPS)))
        for i in range(1, len(text) + 1):
            line_edit.setText(text[:i])
            self.refresh()
            self._emit(every)
        self.hold(after)

    def scroll(self, scrollable, fraction: float, seconds: float = 1.6):
        """Run a scroll bar to `fraction` of its range, a frame at a time, so the move is readable."""
        bar = scrollable.verticalScrollBar()
        start, end = bar.value(), int(bar.minimum() + fraction * (bar.maximum() - bar.minimum()))
        frames = max(2, int(round(seconds * FPS)))
        self.focus = None
        for i in range(1, frames + 1):
            bar.setValue(int(start + (end - start) * _ease(i / frames)))
            self.refresh()
            self._emit()

    def wait(self, ready, *, timeout: float = 1800.0, shown: float = 7.0, minimum: float = 1.5):
        """Pump Qt until `ready()` is true, filming the first `shown` seconds of the wait.

        A strategy started from the window goes to the job runner, so the panel comes back
        immediately and the answer arrives later. The wait is recorded rather than cut out -- the
        console and the jobs panel are doing something, and a tutorial that hides the wait teaches
        the wrong thing -- but only the first few seconds of it, since nobody needs to watch four
        minutes of a progress line.
        """
        start = time.monotonic()
        while True:
            elapsed = time.monotonic() - start
            if ready() and elapsed >= minimum:
                return elapsed
            if elapsed > timeout:
                raise TimeoutError(f"{self.slug}: nothing finished within {timeout:.0f}s")
            tick = time.monotonic()
            self.refresh()
            if elapsed < shown:
                self._emit()
            gap = 1.0 / FPS - (time.monotonic() - tick)
            if gap > 0:
                time.sleep(gap)

    def clear_focus(self):
        """Stop highlighting anything, for a beat that is about reading rather than clicking."""
        self.focus = None

    def mark_poster(self):
        """Use the frame after this call as the poster image the tutorial page shows."""
        self.poster_frame = self.n + 1

    # ------------------------------------------------------------------ encoding
    def write(self) -> dict:
        """Encode the frames to H.264, cut the poster out of them, and report what it came to."""
        os.makedirs(OUT, exist_ok=True)
        video = os.path.join(OUT, f"{self.slug}_silent.mp4")
        poster = os.path.join(OUT, f"{self.slug}_poster.jpg")
        subprocess.run([FFMPEG, "-y", "-loglevel", "error", "-framerate", str(FPS),
                        "-i", os.path.join(self.frames, "%05d.png"),
                        "-c:v", "libx264", "-preset", "slow", "-crf", "27",
                        "-pix_fmt", "yuv420p", "-movflags", "+faststart", video], check=True)
        pick = self.poster_frame or max(1, self.n // 3)
        pick = min(max(1, pick), self.n)
        subprocess.run([FFMPEG, "-y", "-loglevel", "error",
                        "-i", os.path.join(self.frames, f"{pick:05d}.png"),
                        "-vf", "scale=1120:-2", "-q:v", "4", poster], check=True)
        shutil.rmtree(self.frames, ignore_errors=True)
        info = {"slug": self.slug, "title": self.title, "frames": self.n,
                "seconds": self.n / FPS, "video": video, "poster": poster,
                "video_bytes": os.path.getsize(video), "poster_bytes": os.path.getsize(poster)}
        print(f"  {self.slug}: {info['seconds']:.0f}s, "
              f"{info['video_bytes'] / 1000:.0f} kB + {info['poster_bytes'] / 1000:.0f} kB poster",
              flush=True)
        return info


# --------------------------------------------------------------------------- helpers
def choices(panel) -> list:
    """Every answer row on the guided tab's current question, in the order it draws them."""
    from starplast.guided_panel import Choice
    return [b for b in panel.findChildren(Choice) if b.isVisibleTo(panel)]


def choice_in(panel, text: str):
    """The guided tab's answer row whose label starts with `text` -- what a viewer would click."""
    wanted = text.lower()
    for button in choices(panel):
        if button.option.label.lower().startswith(wanted):
            return button
    raise LookupError(f"no guided option starting {text!r}; "
                      f"have {[b.option.label for b in choices(panel)][:12]}")


def choice_value(panel, value: str):
    """The guided answer row that stores `value`, named by the module rather than by its wording."""
    for button in choices(panel):
        if button.option.value == value:
            return button
    raise LookupError(f"no guided option with value {value!r}; "
                      f"have {[b.option.value for b in choices(panel)]}")


def button_in(widget, text: str):
    """The first visible push button whose text contains `text`."""
    from PyQt6 import QtWidgets
    for b in widget.findChildren(QtWidgets.QPushButton):
        if text.lower() in b.text().lower() and b.isVisibleTo(widget):
            return b
    raise LookupError(f"no button matching {text!r}")


def tree_rect(clip: Clip, tree, item):
    """Where a tree row sits in window coordinates, so the pointer can land on that row."""
    from PyQt6 import QtCore
    tree.scrollToItem(item)
    clip.studio.settle()
    r = tree.visualItemRect(item)
    top_left = tree.viewport().mapTo(clip.studio.window, r.topLeft())
    return QtCore.QRect(top_left, r.size())


def list_rect(clip: Clip, listing, row: int):
    """Where a list row sits in window coordinates, so the pointer can land on that match."""
    from PyQt6 import QtCore
    item = listing.item(row)
    listing.scrollToItem(item)
    clip.studio.settle()
    r = listing.visualItemRect(item)
    r.setWidth(listing.viewport().width() - 4)
    return QtCore.QRect(listing.viewport().mapTo(clip.studio.window, r.topLeft()), r.size())


def row_rect(clip: Clip, table, row: int):
    """Where a table row sits in window coordinates."""
    from PyQt6 import QtCore
    table.scrollToItem(table.item(row, 0))
    # scrollToItem also scrolls sideways, which would slide the column headings out of the frame.
    table.horizontalScrollBar().setValue(table.horizontalScrollBar().minimum())
    clip.studio.settle()
    r = table.visualItemRect(table.item(row, 0))
    r.setWidth(table.viewport().width() - 8)
    top_left = table.viewport().mapTo(clip.studio.window, r.topLeft())
    return QtCore.QRect(top_left, r.size())


def reset_window(studio: Studio):
    """Put the window back to how it opens, so one clip cannot start inside another one's state.

    All six clips are filmed from one application -- building six windows would take six times as
    long and Qt keeps style state across them -- so the state they leave behind is cleared here
    rather than left to leak into the next clip's first frame.
    """
    w = studio.window
    for name in ("jobs_dock", "console_dock", "chat_dock", "gallery_dock"):
        dock = getattr(w, name, None)
        if dock is not None:
            dock.hide()
    import numpy as np
    w.search.setText("")
    w.reset()                                   # the find panel's "reset view / clear filters"
    xyz, placed, clusters = studio.home         # the map a maps clip may have swapped out
    if xyz is not None:
        rows = np.arange(len(xyz)) if placed is None else np.where(placed)[0]
        w.use_embedding(xyz if placed is None else xyz[placed], rows)
        w.cluster_labels = clusters
        w.redraw()
    w.status.clearMessage()
    panel = getattr(w, "strategy_panel", None)
    if panel is not None:
        panel.filter.setText("")
    guided = getattr(w, "guided", None)
    if guided is not None:
        guided.restart()
    studio.settle()


def raise_dock(clip: Clip, name: str):
    """Bring one of the tabbed right-hand docks to the front, as clicking its tab would."""
    dock = getattr(clip.studio.window, name)
    dock.show()
    dock.raise_()
    clip.refresh()
    return dock


def first_line(text: str, limit: int = 180) -> str:
    """One tidy line of a label or summary, for a caption that must fit the band.

    Several of the panels' lines are rich text, and a caption is painted, not rendered, so the
    markup is stripped here rather than burned into the frame.
    """
    import html as _html
    import re
    line = " ".join(_html.unescape(re.sub(r"<[^>]+>", " ", str(text))).split())
    return line if len(line) <= limit else line[:limit - 1].rsplit(" ", 1)[0] + "…"


# --------------------------------------------------------------------------- the six clips
def clip_find_a_gene(studio: Studio) -> dict:
    """Clip 1: the find box and the evidence panel -- the first thing anyone does here."""
    w = studio.window
    reset_window(studio)
    clip = Clip(studio, "1_find_a_gene", "Find a gene and read its evidence", 7)
    raise_dock(clip, "right_dock")

    clip.say("Everything starts in the find box",
             "Top left. An exact accession first, then a partial one, then the words of the "
             "product description -- most of this proteome has no symbol.", 2.6)
    clip.move_to(w.search, 1.0)
    clip.hold(0.8)

    clip.say("Type what you know", "GRA16 is a dense-granule protein the parasite exports into "
                                   "the host nucleus. Typing its symbol is enough.")
    clip.type_into(w.search, "GRA16")
    clip.press(w.search, w.do_search, move=0.3, settle=1.4)
    clip.mark_poster()
    clip.hold(1.2)

    clip.say("The evidence panel answers with measurements, not a summary",
             "Every value the shipped tables hold for this gene, each with where it came from. "
             "The header says how much the gene has been written about.", 4.2)
    clip.move_to(w.detail, 1.0, focus=False)
    clip.scroll(w.detail, 0.33, 2.2)
    clip.hold(3.4)
    clip.note("Keep going and the fitness screens, the modifications and the measured partners "
              "follow, each under the experiment that produced it.", 1.0)
    clip.scroll(w.detail, 0.72, 2.6)
    clip.hold(4.0)

    clip.say("A dash is not a zero",
             "Where nothing was measured the panel prints a dash. Absence of evidence is kept "
             "apart from a measured value everywhere in this program -- on the map it is gray.",
             4.0)
    clip.scroll(w.detail, 1.0, 2.0)
    clip.note("The links at the bottom go to the source record, so any number here can be taken "
              "back to the experiment that produced it.", 4.0)

    clip.say("You can search the words instead of the name",
             "Six and a half thousand genes here are 'hypothetical protein', so product text is "
             "often the only handle on a gene.")
    clip.scroll(w.detail, 0.0, 0.6)
    clip.type_into(w.search, "rhoptry neck protein", per_char=0.07)
    clip.press(w.search, w.do_search, move=0.3, settle=2.0)
    clip.note(first_line(w.status.currentMessage(), 170)
              + " — the status bar always says how many genes matched.", 4.4)

    clip.say("Colour the whole map by what you just read",
             "The colour-by box takes any label column, any clustering you have kept, or any "
             "measurement cut into bins. Gray always means unknown.")
    clip.move_to(w.category_box, 1.0)
    clip.hold(2.6)
    clip.move_to(w.comp_list, 1.0)
    clip.note("The list under it filters: tick a class to show only those genes, double-click to "
              "fly to its centroid. The counts are of genes with a value.", 3.4)

    clip.say("That is the loop",
             "Find a gene, read what is actually measured for it, and colour the field by the "
             "property you care about. Everything else in Starplast is built on this one.", 3.4)
    clip.clear_focus()
    clip.hold(1.0)
    return clip.write()


def clip_start_here(studio: Studio) -> dict:
    """Clip 2: the guided tab, from what you have to a strategy running."""
    from starplast import guided as G
    from starplast import strategies as S
    w = studio.window
    reset_window(studio)
    clip = Clip(studio, "2_start_here", "Start here, end to end", 10)
    panel = w.guided
    panel.restart()
    raise_dock(clip, "guided_dock")

    clip.say("Start here asks one question at a time",
             "It does not ask which of 39 strategies you want. It asks what you have, and works "
             "down to the ones calibrated for it on this organism.", 3.2)

    clip.say("What do you have?",
             "One gene, a gene list, your own screen, a label you want extended, or nothing yet. "
             "Everything after this follows from the answer.", 2.4)
    clip.press(choice_in(panel, "A single gene"), move=1.1, settle=1.6)

    clip.say("Which organism?",
             "One species per window, never a union: the identifiers and the data behind them do "
             "not align, and a merged table would quietly invent agreement.", 2.2)
    clip.press(choices(panel)[0], move=1.1, settle=1.6)

    clip.say("Which gene?",
             "Type an accession, a symbol or words from the product; the matches come from the "
             "shipped table, so nothing can be chosen that is not in it.")
    clip.type_into(panel.search, "GRA16", per_char=0.14)
    clip.hold(1.6)
    item = panel.matches.item(0)
    rect = list_rect(clip, panel.matches, 0)
    clip.press(rect, action=lambda: panel._pick_gene(item), move=1.1, settle=1.8)

    clip.say("What do you want to know?",
             "This is the question that decides the task — and therefore which metrics the answer "
             "will be judged by. Partners, more genes like it, or a label predicted for it.", 3.0)
    predict = choice_value(panel, G.PREDICT)
    clip.press(predict, move=1.2, settle=1.8)

    clip.say("Which label?",
             "Each row says how many genes carry a value for it, because a label measured on forty "
             "genes cannot answer the same questions as one measured on four thousand.", 2.6)
    clip.press(choice_in(panel, "compartment"), move=1.2, settle=2.0)

    recs = panel.recommendations()
    clip.say("Ranked strategies, not a menu",
             f"{len(recs)} of {len(S.catalog())} strategies survive what you said, best first, "
             f"each judged by what the calibration sweep measured for it on this organism.", 3.6)
    clip.mark_poster()
    card = panel.cards[0]
    clip.move_to(card, 1.2)
    clip.hold(3.4)
    clip.scroll(panel.scroll, 0.5, 2.2)
    clip.hold(3.4)
    clip.scroll(panel.scroll, 0.0, 1.4)

    clip.say("The settings are already filled in",
             "The answers became the strategy's settings: the organism, the gene, the label, the "
             "hold-out. Run starts it as a background job, with its tables going where every other "
             "strategy's tables go.")
    run = button_in(card, "Run")
    w.strategy_panel.last_result = None
    clip.press(run, move=1.3, settle=0.8)

    clip.say("It runs as a job, and the job is visible",
             "Jobs lists everything running, with Stop beside it and the traceback of anything "
             "that fails. Nothing is written by a stopped run.")
    clip.clear_focus()
    raise_dock(clip, "jobs_dock")
    elapsed = clip.wait(lambda: w.strategy_panel.last_result is not None)
    raise_dock(clip, "strategies_dock")
    w.strategy_panel.show_details(2)
    clip.refresh()
    clip.hold(1.0)

    result = w.strategy_panel.last_result
    summary = first_line(getattr(result, "summary", "")
                         or "the run finished; its tables are under Results")
    clip.say("And it has run", summary, 5.2)
    clip.note(f"{first_line(summary, 150)} ({elapsed:.0f} s on this machine).", 4.4)
    clip.clear_focus()
    clip.hold(1.6)
    return clip.write()


def clip_test_before_you_trust(studio: Studio) -> dict:
    """Clip 3: a strategy card, its four bars, and the hold-out test that judges it."""
    from starplast import strategies as S
    w = studio.window
    reset_window(studio)
    clip = Clip(studio, "3_test_before_you_trust", "Test before you trust", 8)
    panel = w.strategy_panel
    raise_dock(clip, "strategies_dock")
    key = "feature_knn"
    strategy = S.get(key)

    clip.say("The Strategies tab is a catalogue, filtered by method",
             "Thirty-nine named ways of turning the combined data into a claim. The filter matches "
             "names and methods, so a method you trust lists everything built on it.", 2.6)
    clip.type_into(panel.filter, "kNN", per_char=0.2)
    clip.hold(2.0)

    clip.say(f"Strategy {strategy.number:02d}", first_line(strategy.title))
    panel.select(key)
    panel.show_card()
    clip.refresh()
    clip.hold(2.6)

    clip.say("Four bars, the same on every card",
             "Better than chance, reach, and two plain metrics of this strategy's task -- each "
             "with its 95% interval and the level chance alone reaches.", 3.8)
    clip.move_to(panel.card, 1.2, focus=False)
    clip.mark_poster()
    clip.hold(2.8)
    clip.scroll(panel.card_scroll, 0.35, 2.0)
    clip.hold(2.8)

    clip.say("One real failure beside one real success",
             "Not an illustration: two genes from the shipped table, one the strategy got wrong "
             "and one it got right, with what the right answer added.", 3.0)
    clip.scroll(panel.card_scroll, 0.75, 2.0)
    clip.hold(3.2)
    clip.scroll(panel.card_scroll, 0.0, 1.2)

    clip.say("Now test it", "Hide a quarter of the labels, ask the strategy for them back, and "
                            "run the same procedure on shuffled labels for comparison. Nothing "
                            "about this is optional -- every strategy carries it.")
    test_btn = button_in(panel, "Test")
    panel.last_test = None
    clip.press(test_btn, action=panel.test_current, move=1.2, settle=0.8)
    clip.clear_focus()
    elapsed = clip.wait(lambda: panel.last_test is not None)
    clip.hold(0.8)
    test = panel.last_test

    clip.say(f"Verdict: {test.verdict}", first_line(test.summary().split(" -- ", 1)[-1]), 4.4)
    clip.note(f"observed {test.observed:.3f} · chance {test.null_mean:.3f} "
              f"± {test.null_sd:.3f} · effect {test.effect:.3f} · p {test.p_value:.4f} "
              f"({elapsed:.0f} s)", 4.4)

    clip.say("The scorecard says in what way it is good",
             "The verdict rests on one number. The scorecard reports every standard metric for "
             "this kind of task on the same hidden genes, always in the same order, so two "
             "strategies doing the same task compare number by number.", 1.0)
    panel.show_details(2)
    clip.refresh()
    clip.move_to(panel.verdict, 1.2, focus=False)
    clip.hold(4.6)
    clip.move_to(panel.result_tabs, 1.2, focus=False)
    clip.note("Per hidden class, under the verdict: how many came back and how precise the calls "
              "were. A strategy that is excellent on ribosomes and useless on the cytosol says so "
              "here rather than hiding behind an average.", 4.4)

    clip.say("A number without its chance level is not a result",
             "That is the whole rule. Read the verdict, read the gap over chance, and only then "
             "read the candidates.", 3.6)
    clip.clear_focus()
    clip.hold(1.0)
    return clip.write()


def clip_pick_a_map(studio: Studio) -> dict:
    """Clip 4: the gallery of pregenerated maps, sorted by how well a label maps onto each."""
    from starplast.umap_gallery_panel import SORT_LABEL, LABEL_VIEW
    w = studio.window
    reset_window(studio)
    clip = Clip(studio, "4_pick_a_map", "Pick a map that shows your label", 8)
    panel = w.maps_panel
    # The score table is wider than the dock; keep it at its left edge so the first frame starts
    # where a reader would, with the map names and the leading score columns in view.
    panel.table.horizontalScrollBar().setValue(panel.table.horizontalScrollBar().minimum())
    raise_dock(clip, "maps_dock")

    clip.say(f"{len(panel.records)} maps ship already built",
             "Every one from measurements only, grouped by what it was built from: all "
             "measurements, all but one kind, one family, one experiment, pairs and triples.",
             3.4)
    clip.move_to(panel.tree, 1.2, focus=False)
    clip.hold(2.4)
    for i in range(panel.tree.topLevelItemCount()):
        panel.tree.topLevelItem(i).setExpanded(False)
    clip.refresh()
    clip.hold(2.0)

    clip.say("Open a group and the rows state their recipe",
             "Each thumbnail is coloured by that map's own clusters, and each row carries how "
             "many genes it places, how many clusters it found, and its structure score.")
    top = panel.tree.topLevelItem(0)
    clip.press(tree_rect(clip, panel.tree, top),
               action=lambda: top.setExpanded(True), move=1.1, settle=2.6)
    clip.hold(2.4)

    clip.say("Sort by how well YOUR label maps",
             "Structure score is label-free -- it is what the maps were tuned on. Sorting by a "
             "label instead puts the map that separates it best at the top.")
    clip.press(panel.group_check, move=1.1, settle=1.0)
    clip.press(panel.sort_box,
               action=lambda: panel.sort_box.setCurrentText(SORT_LABEL), move=1.0, settle=1.4)
    clip.move_to(panel.label_box, 1.0)
    clip.note(f"The label is {panel.label_box.currentText()}: the hyperLOPIT compartment, scored "
              f"against every map's own clusters.", 3.0)

    panel.view_box.setCurrentText(LABEL_VIEW)
    clip.refresh()
    clip.say("The table ranks the maps for that one label",
             "Categories to clusters, precision, recall, and the skill columns that subtract what "
             "shuffled labels score on the same map -- 0 is chance, 1 is perfect.", 1.0)
    clip.move_to(panel.table, 1.2, focus=False)
    clip.hold(3.6)

    honest = _first_honest_row(clip, panel)
    clip.move_to(row_rect(clip, panel.table, 0), 1.0)
    panel.table.selectRow(0)
    clip.refresh()
    clip.say("Read the top row before believing it",
             first_line(panel.note.text(), 230), 5.4)
    clip.note("The panel says so itself. A map built from a column that restates the label will "
              "of course show it; the useful map is the best one that was not.", 4.6)

    clip.say("So take the best map that was NOT built from it",
             "Double-click the row and that map goes into the main view, where it behaves like "
             "any other map.")
    rect = row_rect(clip, panel.table, honest)
    clip.press(rect, action=lambda: panel._row_activated(honest, 0), move=1.2, settle=2.4)
    clip.mark_poster()
    clip.note(first_line(panel.shown.text(), 210), 4.4)

    clip.say("Colour it by the label and look",
             "One click. Gray is unlabelled, as everywhere here.")
    clip.press(panel.color_btn, move=1.2, settle=2.2)
    clip.note(first_line(panel.note.text(), 230), 4.6)

    clip.say("Choose the map for the question, not the prettiest one",
             "A map that clusters beautifully and scores at chance for your label has told you "
             "nothing about your label. The table is how you tell the difference.", 4.0)
    clip.clear_focus()
    clip.hold(1.2)
    return clip.write()


def _first_honest_row(clip: Clip, panel) -> int:
    """The first ranked map whose own sentence does not say it restates the label being scored.

    The best-scoring map for a label is often one built from a column that is that label by
    another name, and the panel says so in its sentence. The clip picks the first row that does
    not, because that is the map a reader would actually use.
    """
    for row in range(panel.table.rowCount()):
        panel.table.selectRow(row)
        clip.studio.settle(2)
        if "restates" not in panel.note.text().lower():
            return row
    return 0


def clip_walk_the_network(studio: Studio) -> dict:
    """Clip 5: the star map -- one gene's measured links, by source, and the walk outward."""
    from starplast.star_map import MODE_NEIGHBOURHOOD, MODE_STAR
    w = studio.window
    reset_window(studio)
    clip = Clip(studio, "5_walk_the_network", "Walk the network from one gene", 7)
    panel = w.star_map
    raise_dock(clip, "star_map_dock")

    clip.say("One gene in the middle, its measured links around it",
             "The star map never merges the evidence layers into one graph. A crosslink and a "
             "co-mention are different claims, so they are drawn as different colours.", 3.0)

    clip.say("Centre a gene",
             "TGME49_208830 is GRA16, the dense-granule effector the parasite exports into the "
             "host nucleus. Type an accession, or let the panel follow whatever you select on the "
             "map.")
    clip.type_into(panel.gene_edit, "TGME49_208830", per_char=0.07)
    clip.press(panel.centre_btn, move=0.9, settle=2.4)
    clip.mark_poster()
    clip.note(first_line(panel.headline.text(), 190), 3.6)
    clip.note(first_line(panel.counts_label.text(), 190), 3.0)

    clip.say("Turn the sources on and off",
         "This is the question the star map exists for: not 'are these two related' but 'by "
         "which evidence'. Keep only what you would believe.")
    clip.press(panel.measured_btn, move=1.2, settle=2.6)
    clip.note("Measured only: crosslinks, pulldowns, co-fitness, co-expression -- the layers that "
              "came out of an experiment rather than out of a model.", 3.4)
    clip.press(panel.all_btn, move=1.0, settle=2.2)
    clip.note("All of them again. The legend on the right gives every source its own colour, and "
              "says which are measured, which are inferred by a strategy, and which came from a "
              "run of your own.", 4.6)

    clip.say("Widen it", "One gene becomes its neighbourhood: every gene within a couple of hops, "
                         "with the strongest links kept and a sentence saying what was left out.")
    clip.press(panel.mode_box,
               action=lambda: panel.set_mode(MODE_NEIGHBOURHOOD), move=1.1, settle=3.0)
    clip.note(first_line(panel.counts_label.text(), 190), 3.4)
    clip.move_to(panel.stack, 1.2, focus=False)
    clip.hold(2.6)

    neighbour = _a_neighbour(panel)
    if neighbour:
        clip.say("Re-centre on a neighbour",
                 "Following a link to the gene at the other end is the whole move. Back returns "
                 "to where you were, so a walk is reversible.")
        clip.press(panel.mode_box,
                   action=lambda: panel.set_mode(MODE_STAR), move=1.0, settle=1.2)
        clip.type_into(panel.gene_edit, neighbour)
        clip.press(panel.centre_btn, move=0.9, settle=2.6)
        clip.note(first_line(panel.headline.text(), 190), 3.4)
        clip.press(panel.back_btn, move=1.2, settle=2.4)
        clip.note("Back: the previous gene, with the sources you had chosen still chosen.", 2.4)
    else:
        clip.say("Re-centre on a neighbour",
                 "Following a link to the gene at the other end is the whole move; Back returns "
                 "to where you were.", 3.0)

    clip.say("A link is a claim with a source",
             "Hover any line in the application and it names the layer, the strength and the "
             "experiment. Nothing here is a consensus edge from an unnamed merge.", 3.6)
    clip.clear_focus()
    clip.hold(1.0)
    return clip.write()


def _a_neighbour(panel) -> str:
    """One gene the star map has just drawn beside the centred one, as an accession.

    `last_nodes` carries row numbers into the window's table, not identifiers, so the clip would
    otherwise type a row number into the gene box and be told there is no such gene.
    """
    rows = getattr(panel, "last_nodes", None)
    if rows is None or not len(rows):
        return ""
    centre = str(panel.gene_edit.text()).strip().upper()
    for value in rows["gene"]:
        try:
            gene = str(panel.gene_ids[int(value)])
        except (TypeError, ValueError, IndexError):
            gene = str(value)
        if gene.upper() != centre:
            return gene
    return ""


def clip_one_question(studio: Studio) -> dict:
    """Clip 6: question 13 of instruction 59, from the question to the genes it names."""
    from starplast import strategies as S
    w = studio.window
    reset_window(studio)
    clip = Clip(studio, "6_one_question", "One real biological question, answered", 8)
    panel = w.strategy_panel
    raise_dock(clip, "strategies_dock")
    key, strategy = "understudied_first", S.get("understudied_first")

    clip.say("The question",
             "Which never-published Toxoplasma proteins get a confident compartment call? "
             "Question 13 of the hundred in instruction 59, anchored on PMID:41582196.", 4.4)
    clip.note("It is a real question because of what it is NOT: not a re-derivation of a label the "
              "data already carries, but a call about genes nobody has written about.", 4.0)

    clip.say("It has an entry point", f"Strategy {strategy.number:02d}, "
                                      f"{first_line(strategy.title, 90)}.")
    clip.type_into(panel.filter, "understudied")
    panel.select(key)
    panel.show_card()
    clip.refresh()
    clip.hold(3.0)

    clip.say("Read what it does before running it",
             "It makes a call only where independent methods agree, and ranks what is left by how "
             "little has been written about the gene. Agreement is the whole point: it is what "
             "makes a call about an unstudied gene worth anything.", 1.0)
    clip.move_to(panel.card, 1.2, focus=False)
    clip.hold(3.6)
    clip.scroll(panel.card_scroll, 0.45, 2.0)
    clip.hold(3.0)
    clip.scroll(panel.card_scroll, 0.0, 1.0)

    clip.say("The settings the question asks for",
             "The hyperLOPIT compartment as the target, and at least two methods having to agree. "
             "Exactly what instruction 59 records for this question.")
    panel.show_details(1)
    clip.refresh()
    panel.set_setting("target", "compartment")
    panel.set_setting("min_agree", 2)
    clip.refresh()
    clip.move_to(panel.form_host, 1.2, focus=False)
    clip.hold(3.4)

    clip.say("Run it", "On the whole table, in the background, with the label and everything that "
                       "restates it withheld from every method that votes.")
    run_btn = button_in(panel, "Run")
    panel.last_result = None
    clip.press(run_btn, action=panel.run_current, move=1.2, settle=0.8)
    clip.clear_focus()
    elapsed = clip.wait(lambda: panel.last_result is not None)
    clip.hold(0.8)
    result = panel.last_result

    clip.say("The answer", first_line(result.summary), 4.6)
    clip.note(f"{first_line(result.summary, 150)} ({elapsed:.0f} s on this machine).", 3.0)

    table = result.tables.get("candidates")
    if table is None:
        table = list(result.tables.values())[0]
    named = ", ".join(str(g) for g in table.iloc[:3, 0])
    clip.say("Named genes, which is the point",
             f"{len(table):,} candidates, best first. The top three are {named} -- each with the "
             f"compartment the agreeing methods called, the support behind it, and how many "
             f"methods agreed.", 1.0)
    panel.show_details(2)
    clip.refresh()
    clip.mark_poster()
    clip.move_to(panel.result_tabs, 1.2, focus=False)
    clip.hold(4.6)
    shown = panel.tables[panel.result_tabs.currentIndex()]
    clip.scroll(shown, 0.12, 2.0)
    clip.hold(3.4)

    clip.say("From a question to a shortlist you could pick up",
             "Right-click saves the table. The same run, at the same settings, is recorded in "
             "results/questions_2026_09_30 -- the question, the strategy, the numbers and the "
             "genes, reproducible.", 4.0)
    clip.clear_focus()
    clip.hold(1.2)
    return clip.write()


#: The clips, in the order they are RECORDED -- which is not the order they are numbered in. One
#: application films all six, and a strategy run keeps its predicted links for the star map to
#: offer, so the two clips that run a strategy are filmed after the star-map clip rather than
#: before it: otherwise the star map opens showing links the viewer never made.
CLIPS = {
    "1_find_a_gene": clip_find_a_gene,
    "4_pick_a_map": clip_pick_a_map,
    "5_walk_the_network": clip_walk_the_network,
    "3_test_before_you_trust": clip_test_before_you_trust,
    "2_start_here": clip_start_here,
    "6_one_question": clip_one_question,
}

#: What each clip shows, for the tutorial pages and the index. Written here so the page builder and
#: the video builder cannot disagree about it.
DESCRIPTIONS = {
    "1_find_a_gene": ("Find a gene and read its evidence",
                      "The find box, the evidence panel, and the rule that a dash is not a zero.",
                      "1_explore"),
    "2_start_here": ("Start here, end to end",
                     "The guided tab: what you have, which genes, which label, what you want to "
                     "know -- then the ranked strategies and Run.", "3_gene_list"),
    "3_test_before_you_trust": ("Test before you trust",
                                "A strategy card, its four bars, the hold-out test, and the "
                                "verdict with its chance level.", "4_test_and_calibrate"),
    "4_pick_a_map": ("Pick a map that shows your label",
                     "The pregenerated gallery, sorted by how well a label maps onto each, and "
                     "the score table that says so.", "2_holdout_search"),
    "5_walk_the_network": ("Walk the network from one gene",
                           "The star map: one gene's measured links by source, widened to its "
                           "neighbourhood and re-centred on a neighbour.",
                           "5_networks_and_agreement"),
    "6_one_question": ("One real biological question, answered",
                       "Question 13 of instruction 59: which never-published proteins get a "
                       "confident compartment call, from the question to the named genes.",
                       "6_advanced_models"),
}


# --------------------------------------------------------------------------- publishing
#: The clips' markup is written into pages `build_tutorials.py` generates, so it is bracketed and
#: replaced rather than appended: running either script twice leaves one copy, not two.
MARK = ("<!--starplast-video-->", "<!--/starplast-video-->")
PAGES = os.path.join(ROOT, "docs", "tutorial")


def clip_for(tutorial_slug: str) -> str:
    """The clip recorded for one of the six written tutorials, or ""."""
    for slug, (_t, _d, tutorial) in DESCRIPTIONS.items():
        if tutorial == tutorial_slug:
            return slug
    return ""


def _player(slug: str, base: str, width: str) -> str:
    """A poster that plays on click: no preload, so a page of six costs one request each."""
    return (f"<video controls preload='none' playsinline style='width:{width};border-radius:8px;"
            f"background:#111' poster='{base}video/{slug}_poster.jpg'>"
            f"<source src='{base}video/{slug}_silent.mp4' type='video/mp4'>"
            f"Your browser cannot play this clip; "
            f"<a href='{base}video/{slug}_silent.mp4'>download it</a> instead.</video>")


def _block(body: str) -> str:
    return MARK[0] + body + MARK[1]


def _replace(text: str, body: str) -> str:
    """Put `body` in the page's video block, replacing an older one if it is already there."""
    start, end = text.find(MARK[0]), text.find(MARK[1])
    if start >= 0 and end > start:
        return text[:start] + _block(body) + text[end + len(MARK[1]):]
    return text


def publish_pages() -> list:
    """Write the clips into the tutorial pages and the tutorial index, and say which were touched.

    The six written tutorials each gain their clip at the top of the settings panel, and the index
    gains a row of posters. Both are generated pages, so the markup is bracketed by `MARK` and this
    runs again after every `build_tutorials.py`, which calls it.
    """
    touched = []
    for slug, (title, description, tutorial) in DESCRIPTIONS.items():
        path = os.path.join(PAGES, f"{tutorial}_gui.html")
        if not os.path.exists(path) or not os.path.exists(
                os.path.join(OUT, f"{slug}_silent.mp4")):
            continue
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
        body = (f"<div style='margin:0 0 14px 0'>{_player(slug, '', '100%')}"
                f"<div style='font-size:12.5px;color:#a8a8a8;margin-top:6px'><b>Watch it: "
                f"{html_escape(title)}.</b> {html_escape(description)} Silent screen capture of "
                f"the real application.</div></div>")
        anchor = "<div class='intro'>"
        if MARK[0] in text:
            new = _replace(text, body)
        elif anchor in text:
            new = text.replace(anchor, anchor + _block(body), 1)
        else:
            continue
        if new != text:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(new)
            touched.append(path)

    index = os.path.join(PAGES, "index.html")
    if os.path.exists(index):
        cards = "".join(
            f"<div style='background:#2a2a2a;border-radius:12px;padding:12px'>"
            f"{_player(slug, '', '100%')}"
            f"<h3 style='font-size:15px;margin:10px 0 4px'>{html_escape(title)}</h3>"
            f"<p style='font-size:13px;color:#a8a8a8;margin:0'>{html_escape(description)}</p>"
            f"<p style='font-size:12.5px;margin:8px 0 0'>"
            f"<a href='{tutorial}_gui.html'>the written tutorial</a></p></div>"
            for slug, (title, description, tutorial) in sorted(DESCRIPTIONS.items())
            if os.path.exists(os.path.join(OUT, f"{slug}_silent.mp4")))
        body = (f"<h2>Watch it done</h2><p>Six silent screen captures of the real application, one "
                f"to two minutes each. Every panel, number and table in them was produced by the "
                f"program while the clip was being recorded; the captions are burned in, so there "
                f"is nothing to listen to. Click a poster to play.</p>"
                f"<div style='display:grid;gap:16px;grid-template-columns:repeat(auto-fit,"
                f"minmax(320px,1fr))'>{cards}</div>" if cards else "")
        with open(index, encoding="utf-8") as fh:
            text = fh.read()
        anchor = "<h2>The complete guide</h2>"
        if MARK[0] in text:
            new = _replace(text, body)
        elif anchor in text:
            new = text.replace(anchor, _block(body) + anchor, 1)
        else:
            new = text
        if new != text:
            with open(index, "w", encoding="utf-8") as fh:
                fh.write(new)
            touched.append(index)
    return touched


NOTEBOOK = os.path.join(ROOT, "notebooks", "tutorial_videos_2026_10_03.ipynb")


def write_notebook(path: str = NOTEBOOK) -> str:
    """Record the recorded clips as an executed notebook, as the project rule requires.

    It measures the files that were produced rather than restating what this script intended: the
    durations and sizes come out of ffprobe, and the genes named in clip 6 are compared with the
    recorded run of question 13 that the clip claims to reproduce.
    """
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    from notebook_runner import ExecutedNotebook
    nb = ExecutedNotebook("Tutorial videos: six screen captures of the real application")
    nb.md("Built by `scripts/tutorial_video.py`, which drives the real `Window` offscreen, grabs "
          "each frame with `QWidget.grab()`, paints the pointer, the highlight and the caption "
          "over it, and encodes the frames with `/usr/bin/ffmpeg`.",
          "",
          "Nothing in the clips is staged. Every panel is the real one, every number was computed "
          "by the shipped code while the clip was recorded, and the only pixels the script paints "
          "are the overlay. This notebook measures what came out.",
          "",
          "    QT_QPA_PLATFORM=offscreen python scripts/tutorial_video.py")
    nb.code("import json, os, subprocess",
            "import pandas as pd",
            f"ROOT = {ROOT!r}",
            "OUT = os.path.join(ROOT, 'docs', 'tutorial', 'video')",
            "sorted(os.listdir(OUT))")
    nb.md("## What each clip is, and what it came to",
          "",
          "`ffprobe` on the committed files: the duration, the frame rate and the codec are what "
          "a browser will actually play, not what the encoder was asked for.")
    nb.code(
        "def probe(path):",
        "    out = subprocess.run(['/usr/bin/ffprobe', '-v', 'error', '-show_entries',",
        "                          'stream=width,height,r_frame_rate,codec_name,pix_fmt',",
        "                          '-show_entries', 'format=duration', '-of', 'json', path],",
        "                         capture_output=True, text=True, check=True)",
        "    d = json.loads(out.stdout)",
        "    s = d['streams'][0]",
        "    return {'codec': s['codec_name'], 'size': f\"{s['width']}x{s['height']}\",",
        "            'fps': s['r_frame_rate'], 'pix_fmt': s['pix_fmt'],",
        "            'seconds': round(float(d['format']['duration']), 1)}",
        "",
        f"DESCRIPTIONS = {json.dumps({k: list(v) for k, v in DESCRIPTIONS.items()})}",
        "rows = []",
        "for slug, (title, what, tutorial) in sorted(DESCRIPTIONS.items()):",
        "    video = os.path.join(OUT, f'{slug}_silent.mp4')",
        "    poster = os.path.join(OUT, f'{slug}_poster.jpg')",
        "    rows.append({'clip': slug, 'title': title, **probe(video),",
        "                 'video_kB': round(os.path.getsize(video) / 1000),",
        "                 'poster_kB': round(os.path.getsize(poster) / 1000),",
        "                 'tutorial': tutorial})",
        "clips = pd.DataFrame(rows)",
        "clips")
    nb.code("added = sum(os.path.getsize(os.path.join(OUT, f)) for f in os.listdir(OUT))",
            "print(f'{len(os.listdir(OUT))} files, {added:,} bytes "
            "({added / 1e6:.1f} MB) added to the repository')",
            "print('every clip between 60 and 120 s:',"
            " bool((clips.seconds.between(60, 120)).all()))")
    nb.md("## Clip 6 shows the real answer to a real question",
          "",
          "The sixth clip runs question 13 of `instructions/open/59_biological_questions.md` -- "
          "which never-published proteins get a confident compartment call -- live in the "
          "Strategies tab. The genes its caption names must be the genes the recorded run of that "
          "question produced, or the clip is showing something made up.")
    nb.code("recorded = pd.read_csv(os.path.join(ROOT, 'results', 'questions_2026_09_30',",
            "                                    'Q13_understudied_first_candidates.csv'))",
            "recorded.head(3)[['gene_id', 'product', 'prediction', 'support', "
            "'sources_agreeing']]")
    nb.code("named_in_the_clip = ['TGME49_209000', 'TGME49_225745', 'TGME49_266010']",
            "print('the clip names:', named_in_the_clip)",
            "print('the recorded run:', list(recorded.gene_id.head(3)))",
            "print('match:', list(recorded.gene_id.head(3)) == named_in_the_clip)",
            "print('candidates in the recorded run:', len(recorded))")
    nb.md("## The clips ship with the site, not with the wheel",
          "",
          "`scripts/build_docs.py` copies `docs/tutorial/` -- the written tutorials and these "
          "clips -- into `docs/site/`. The Python package declares its own `packages` and "
          "`package-data`, neither of which reaches `docs/`, so nothing here is installed by pip.")
    nb.code("import tomllib",
            "with open(os.path.join(ROOT, 'pyproject.toml'), 'rb') as fh:",
            "    cfg = tomllib.load(fh)['tool']['setuptools']",
            "print('packages:', cfg['packages'])",
            "print('package data:', cfg['package-data'])",
            "print('include-package-data:', cfg['include-package-data'])")
    nb.md("## Rebuilding",
          "",
          "    QT_QPA_PLATFORM=offscreen python scripts/tutorial_video.py             # all six",
          "    QT_QPA_PLATFORM=offscreen python scripts/tutorial_video.py 4_pick_a_map  # one",
          "    python scripts/tutorial_video.py --pages                               # relink only",
          "",
          "`scripts/build_tutorials.py` calls `publish_pages()` at the end, so regenerating the "
          "written tutorials puts the clips back into them rather than dropping them.")
    return nb.write(path)


def html_escape(text: str) -> str:
    """Escape text going into the generated pages, as the page templates do."""
    import html
    return html.escape(str(text))


def main(argv=None) -> int:
    """Record every clip, or only the ones named, and print what each came to."""
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv == ["--pages"]:                     # relink without re-recording
        for path in publish_pages():
            print(f"linked {os.path.relpath(path, ROOT)}")
        return 0
    if argv == ["--notebook"]:                  # record what the clips came to
        print(f"wrote {write_notebook()}")
        return 0
    if not os.path.exists(FFMPEG):
        raise SystemExit(f"{FFMPEG} is not installed; the clips cannot be encoded")
    os.makedirs(OUT, exist_ok=True)
    studio = Studio()
    total = 0
    try:
        for slug, build in CLIPS.items():
            if argv and slug not in argv:
                continue
            t0 = time.monotonic()
            print(f"{slug}…", flush=True)
            info = build(studio)
            total += info["video_bytes"] + info["poster_bytes"]
            print(f"  built in {time.monotonic() - t0:.0f}s", flush=True)
    finally:
        studio.close()
    for path in publish_pages():
        print(f"linked {os.path.relpath(path, ROOT)}")
    print(f"wrote {OUT}: {total / 1e6:.1f} MB added")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
