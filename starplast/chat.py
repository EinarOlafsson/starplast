"""An assistant pane that can see the map.

The same shape as the AI console in spaCR: shell out to whichever vendor coding-agent CLI the user
already has logged in (`claude`, `codex`, `gemini`) rather than asking for an API key. Those CLIs
authenticate against a chat subscription, so there is no key to store, leak, or forget to revoke, and
the one in this repository stays a viewer rather than becoming a thing that bills people.

What makes it worth having here rather than in a browser tab is context. The panel builds a briefing
from the actual state -- the selected gene and its evidence, the coloring, how many genes are
filtered out, and the standing caveats about what this map does and does not recover -- so a question
about "this gene" has a referent. An assistant that cannot see the screen answers about Toxoplasma in
general, which the user can already get anywhere.

The caveats are not decoration. The held-out search says localization is recovered worse than study
effort, so an assistant left to free-associate about a gene's neighbours would produce exactly the
overclaim the rest of the project is built to prevent.
"""
from __future__ import annotations

import queue
import shutil
import subprocess
import threading
import time
from html import escape as html_escape
from typing import Iterator, Optional

from PyQt6 import QtCore, QtWidgets

# The standing caveats, sent with every question. Kept here as one string so there is a single place
# to correct when the numbers change, rather than a claim reconstructed per prompt.
GROUNDING = """\
You are answering inside starplast, a gene-evidence browser for Toxoplasma and Plasmodium.
Use the organism, source measurements and analysis supplied with this question.

The packaged maps use a saved balanced feature recipe with numeric biological measurements,
median imputation and robust scaling. Categorical compartment labels are not one-hot encoded;
numeric localization confidence can still be an input. No literature column is an input, but
measurement coverage and assay selection can still shape the map. The archive records the
actual features, algorithm, backend and ordered gene IDs. An interactive recipe can differ.

The display is exploratory, not an independently validated function prediction. Proximity alone
does not establish localization, interaction or mechanism. Guided predictions use separate
group-held-out evaluation, target-derived input exclusions and calibration where labels permit.
Use the actual exported run for performance claims; do not reuse historical cluster-recovery
figures as if they measured a newly fitted model. Distinguish measured evidence, transferred
annotations and model hypotheses. Unsupported genes receive no call; grey means unknown.

Describe the evidence and its source. Say when a requested inference is unsupported. Prefer
"these genes are near each other in this embedding" over "these genes are related". Be brief."""


class Provider:
    """A vendor CLI that can answer a prompt non-interactively."""

    def __init__(self, name: str, label: str, cli: str, hint: str, login: str):
        self.name, self.label, self.cli = name, label, cli
        self.install_hint, self.login_command = hint, login

    def available(self) -> bool:
        """Whether this provider's CLI is on PATH."""
        return shutil.which(self.cli) is not None

    def argv(self, prompt: str, system: str) -> list[str]:
        """The command line that asks this CLI one question non-interactively."""
        if self.name == "claude":
            return [self.cli, "-p", prompt, "--append-system-prompt", system]
        if self.name == "codex":
            return [self.cli, "exec", f"{system}\n\n{prompt}"]
        return [self.cli, "-p", f"{system}\n\n{prompt}"]


PROVIDERS = [
    Provider("claude", "Claude (via Claude Code)", "claude",
             "curl -fsSL https://claude.ai/install.sh | bash", "claude setup-token"),
    Provider("codex", "ChatGPT (via Codex CLI)", "codex",
             "npm install -g @openai/codex", "codex login"),
    Provider("gemini", "Gemini (via Gemini CLI)", "gemini",
             "npm install -g @google/gemini-cli", "gemini auth login"),
]


def available_providers() -> list[Provider]:
    """Providers whose CLI is actually on PATH, in preference order."""
    return [p for p in PROVIDERS if p.available()]


def stream(provider: Provider, prompt: str, system: str, timeout: int = 180,
           register=None) -> Iterator[str]:
    """Yield stdout chunks from the CLI. Errors are yielded as text, never raised.

    This runs on a worker thread feeding a widget, and an exception there is either a silent dead
    panel or a crash. A message the user can read is strictly better than both.

    Reading happens on a helper thread feeding a queue, rather than by iterating readline() here.
    That looks like the long way round and is the only version where `timeout` means anything: a CLI
    that hangs without printing never returns from readline, so a timeout applied to the wait()
    afterwards is never reached and the panel waits forever on a process that will never speak.

    `register` receives the Popen as soon as it exists, so a caller can kill it. Without that, "stop"
    can only set a flag that is checked between chunks -- and a thread blocked on a silent process
    never reaches the check, so closing the window would leave a live thread and abort the process.
    """
    argv = provider.argv(prompt, system)
    try:
        proc = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                text=True, bufsize=1)
    except FileNotFoundError:
        yield f"`{provider.cli}` is not on PATH.\n\nInstall it with:\n  {provider.install_hint}"
        return
    except Exception as exc:
        yield f"could not start `{provider.cli}`: {type(exc).__name__}: {exc}"
        return
    if register is not None:
        register(proc)

    q: "queue.Queue[str | None]" = queue.Queue()

    def reader():
        try:
            for line in iter(proc.stdout.readline, ""):
                q.put(line)
        except Exception:
            pass                       # the pipe closing under a kill is expected, not an error
        finally:
            q.put(None)

    t = threading.Thread(target=reader, daemon=True)
    t.start()
    deadline = time.monotonic() + timeout
    timed_out = False
    try:
        while True:
            left = deadline - time.monotonic()
            if left <= 0:
                timed_out = True
                proc.kill()
                yield f"\n\n`{provider.cli}` timed out after {timeout}s."
                break
            try:
                item = q.get(timeout=min(left, 0.2))
            except queue.Empty:
                continue
            if item is None:
                break
            yield item
        if not timed_out:
            try:
                proc.wait(timeout=max(deadline - time.monotonic(), 0.1))
            except subprocess.TimeoutExpired:
                proc.kill()
            if proc.returncode:
                err = (proc.stderr.read() or "").strip()
                # A non-zero exit almost always means "not logged in", and the fix is a command.
                yield (f"\n\n`{provider.cli}` exited {proc.returncode}."
                       f"{chr(10) + err if err else ''}"
                       f"\n\nIf this is an authentication failure, run:\n  {provider.login_command}")
    finally:
        if proc.poll() is None:
            proc.kill()
        for s in (proc.stdout, proc.stderr):
            try:
                s.close()
            except Exception:
                pass


class _Worker(QtCore.QThread):
    """Streams one reply on a background thread, and can be stopped mid-stream."""
    chunk = QtCore.pyqtSignal(str)
    done = QtCore.pyqtSignal()

    def __init__(self, provider, prompt, system, parent=None):
        super().__init__(parent)
        self.provider, self.prompt, self.system = provider, prompt, system
        self._stop = False
        self._proc = None
        self._lock = threading.Lock()

    def stop(self):
        """Ask the worker to stop, and kill the child so it actually can.

        The flag alone is checked only between chunks, so a thread blocked on a CLI that has not
        printed anything never reaches it -- and a QThread still running when its panel is destroyed
        aborts the interpreter. Killing the child unblocks the read.
        """
        with self._lock:
            self._stop = True
            proc = self._proc
        self._kill(proc)

    def _register(self, proc):
        """Called from the worker thread the moment the child exists.

        Under the lock, and it kills immediately if a stop already arrived: closing the panel in the
        instant between start() and the process existing would otherwise find `_proc` still None,
        skip the kill, and leave the thread reading a process nobody is waiting for. That race is
        exactly the case that matters -- someone closing the window as a reply begins.
        """
        with self._lock:
            self._proc = proc
            stopped = self._stop
        if stopped:
            self._kill(proc)

    @staticmethod
    def _kill(proc):
        if proc is not None and proc.poll() is None:
            try:
                proc.kill()
            except Exception:
                pass

    def run(self):
        for piece in stream(self.provider, self.prompt, self.system, register=self._register):
            if self._stop:
                break
            self.chunk.emit(piece)
        self.done.emit()


class _Input(QtWidgets.QTextEdit):
    """Enter sends, Shift+Enter makes a newline -- the convention every chat box uses."""

    submitted = QtCore.pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setPlaceholderText("Ask about the selected gene, or what you are looking at…")
        self.setMaximumHeight(84)

    def keyPressEvent(self, ev):
        if (ev.key() in (QtCore.Qt.Key.Key_Return, QtCore.Qt.Key.Key_Enter)
                and not (ev.modifiers() & QtCore.Qt.KeyboardModifier.ShiftModifier)):
            self.submitted.emit()
            return
        super().keyPressEvent(ev)


class ChatPanel(QtWidgets.QWidget):
    """The assistant pane. `context_provider` is a callable returning a string about the map."""

    def __init__(self, context_provider=None, parent=None):
        super().__init__(parent)
        self.context_provider = context_provider or (lambda: "")
        self._worker: Optional[_Worker] = None
        self._reply = ""

        L = QtWidgets.QVBoxLayout(self)
        L.setContentsMargins(6, 6, 6, 6)
        L.setSpacing(6)

        bar = QtWidgets.QHBoxLayout()
        self.provider_box = QtWidgets.QComboBox()
        self.refresh_providers()
        self.provider_box.setToolTip(
            "These shell out to a coding-agent CLI you have already logged in, so there is no API "
            "key to store and nothing here is metered separately.")
        bar.addWidget(self.provider_box, 1)
        self.send_btn = QtWidgets.QPushButton("send")
        self.send_btn.setToolTip(
            "Send the question, along with a briefing on what is currently on screen. Enter sends, "
            "shift+enter starts a new line. The assistant is told not to infer a gene's function "
            "from its neighbours, because this map does not support that.")
        self.send_btn.clicked.connect(self.send)
        bar.addWidget(self.send_btn)
        L.addLayout(bar)

        self.log = QtWidgets.QTextBrowser()
        self.log.setOpenExternalLinks(True)
        self.log.setHtml("<p style='color:#888'>Ask a question. The assistant is told what is on "
                         "screen, and is told not to infer function from position.</p>")
        L.addWidget(self.log, 1)

        self.input = _Input()
        self.input.submitted.connect(self.send)
        L.addWidget(self.input)

    def refresh_providers(self):
        """Repopulate the provider list from what is actually installed."""
        self.provider_box.clear()
        got = available_providers()
        for p in got:
            self.provider_box.addItem(p.label, p.name)
        if not got:
            self.provider_box.addItem("no assistant CLI found", None)
            self.provider_box.setEnabled(False)

    def current_provider(self) -> Optional[Provider]:
        """The selected provider, or None when nothing is installed."""
        name = self.provider_box.currentData()
        return next((p for p in PROVIDERS if p.name == name), None)

    def system_prompt(self) -> str:
        """The standing caveats plus a briefing on what is currently on screen."""
        ctx = ""
        try:
            ctx = self.context_provider() or ""
        except Exception as exc:
            # Context is a nicety; failing to build it must not block the question.
            ctx = f"(context unavailable: {type(exc).__name__})"
        return f"{GROUNDING}\n\nCurrent state of the viewer:\n{ctx}" if ctx else GROUNDING

    def _append(self, who: str, text: str):
        # Escaped, because a model reply containing "<" would otherwise be parsed as markup and
        # silently swallow the rest of the answer.
        safe = html_escape(text).replace("\n", "<br>")
        color = self._accent(who)
        self.log.append(f"<p style='margin:6px 0'><b style='color:{color}'>{html_escape(who)}</b>"
                        f"<br>{safe}</p>")

    def _accent(self, who: str) -> str:
        """Two distinguishable colors from the active theme rather than fixed hexes.

        Hardcoded chat colors are the usual way a panel becomes unreadable on a light theme.
        """
        try:
            from . import theme as TH
            p = TH.palette_for(getattr(self.window(), "theme", "dark"))
            return p["accent"] if who == "you" else p.get("accent_lo", p["fg"])
        except Exception:
            return "#7aa2f7" if who == "you" else "#9ece6a"

    def send(self):
        """Send the question in the input box, streaming the reply."""
        prompt = self.input.toPlainText().strip()
        if not prompt or self._worker is not None:
            return
        provider = self.current_provider()
        if provider is None:
            hints = "\n".join(f"  {p.label}: {p.install_hint}" for p in PROVIDERS)
            self.log.append("<p style='color:#e06c75'>No assistant CLI on PATH. Install one:</p>"
                            f"<pre>{hints}</pre>")
            return
        self.input.clear()
        self._append("you", prompt)
        self._reply = ""
        self.send_btn.setText("stop")
        self.send_btn.clicked.disconnect()
        self.send_btn.clicked.connect(self.stop)

        self._worker = _Worker(provider, prompt, self.system_prompt(), self)
        self._worker.chunk.connect(self._on_chunk)
        self._worker.done.connect(self._on_done)
        self._worker.start()

    def _on_chunk(self, text: str):
        self._reply += text

    def _on_done(self):
        self._append(self.provider_box.currentText().split(" (")[0],
                     self._reply.strip() or "(no output)")
        self._worker = None
        self.send_btn.setText("send")
        self.send_btn.clicked.disconnect()
        self.send_btn.clicked.connect(self.send)

    def stop(self):
        """Stop a reply in progress, killing the child process so it actually stops."""
        if self._worker is not None:
            self._worker.stop()

    def closeEvent(self, ev):
        self.stop()
        if self._worker is not None:
            self._worker.wait(2000)
        super().closeEvent(ev)
