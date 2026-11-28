"""Tests for the --web helpers in main.py."""

import sys
import threading
import unittest
import urllib.request
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import main  # noqa: E402


class TestWebUrl(unittest.TestCase):
    def test_localhost(self):
        self.assertEqual(main.web_url("127.0.0.1", 8000), "http://127.0.0.1:8000/web/")

    def test_wildcard_shown_as_localhost(self):
        self.assertEqual(main.web_url("0.0.0.0", 8080), "http://127.0.0.1:8080/web/")
        self.assertEqual(main.web_url("::", 8080), "http://[::1]:8080/web/")

    def test_specific_address(self):
        self.assertEqual(main.web_url("192.168.1.20", 8000), "http://192.168.1.20:8000/web/")
        self.assertEqual(main.web_url("::1", 8000), "http://[::1]:8000/web/")


class TestHostOption(unittest.TestCase):
    def parse(self, *argv):
        with mock.patch.object(sys, "argv", ["main.py", *argv]), \
             mock.patch.object(main, "serve_web") as serve:
            main.main()
        return serve

    def test_default_host_is_local_only(self):
        self.parse("--web").assert_called_once_with(main.DEFAULT_PORT, "127.0.0.1")

    def test_custom_host_and_port(self):
        self.parse("--web", "--host", "0.0.0.0", "--port", "9000").assert_called_once_with(9000, "0.0.0.0")


class TestServeWeb(unittest.TestCase):
    def test_serves_site_without_opening_browser_on_non_local_host(self):
        # Bind to all interfaces on a free port and fetch the page.
        started = threading.Event()
        servers = []
        real_class = main.http.server.ThreadingHTTPServer

        class Recording(real_class):
            def __init__(self, *a, **kw):
                super().__init__(*a, **kw)
                servers.append(self)
                started.set()

        with mock.patch.object(main.http.server, "ThreadingHTTPServer", Recording), \
             mock.patch.object(main.webbrowser, "open") as open_browser, \
             mock.patch("builtins.print"):
            t = threading.Thread(target=main.serve_web, args=(0, "0.0.0.0"), daemon=True)
            t.start()
            self.assertTrue(started.wait(5))
            port = servers[0].server_address[1]
            try:
                with urllib.request.urlopen(f"http://127.0.0.1:{port}/web/index.html", timeout=5) as r:
                    self.assertEqual(r.status, 200)
            finally:
                servers[0].shutdown()
                t.join(5)
        open_browser.assert_not_called()


if __name__ == "__main__":
    unittest.main()
