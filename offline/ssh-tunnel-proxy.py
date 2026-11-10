#!/usr/bin/env python3
"""SSH ProxyCommand helper: tunnels through the runtime HTTP proxy (port 3130 = Tailscale).
Uses raw fd IO to avoid buffered-stdin shutdown crashes."""
import socket, sys, os, urllib.parse, base64, select

def main():
    dest_host, dest_port = sys.argv[1], int(sys.argv[2])
    u = urllib.parse.urlparse(os.environ["HTTPS_PROXY"])
    proxy_host, proxy_port = u.hostname, 3130
    auth = base64.b64encode(
        f"{urllib.parse.unquote(u.username or '')}:{urllib.parse.unquote(u.password or '')}".encode()
    ).decode()
    s = socket.create_connection((proxy_host, proxy_port), timeout=30)
    s.sendall((
        f"CONNECT {dest_host}:{dest_port} HTTP/1.1\r\n"
        f"Host: {dest_host}:{dest_port}\r\n"
        f"Proxy-Authorization: Basic {auth}\r\n\r\n"
    ).encode())
    buf = b""
    while b"\r\n\r\n" not in buf:
        chunk = s.recv(1024)
        if not chunk:
            sys.stderr.write("proxy closed during CONNECT\n")
            os._exit(1)
        buf += chunk
    status_line = buf.split(b"\r\n")[0].decode(errors="replace")
    if "200" not in status_line:
        sys.stderr.write(f"CONNECT failed: {status_line}\n")
        os._exit(1)
    extra = buf.split(b"\r\n\r\n", 1)[1]
    if extra:
        os.write(1, extra)

    s.setblocking(False)
    fds = [0, s.fileno()]
    bufs = {s.fileno(): b""}
    open_in = True
    while True:
        watch = [s.fileno()] + ([0] if open_in else [])
        r, _, _ = select.select(watch, [], [], 60)
        if not r:
            continue
        for fd in r:
            if fd == 0:
                data = os.read(0, 4096)
                if not data:
                    open_in = False
                    try:
                        s.shutdown(socket.SHUT_WR)
                    except Exception:
                        pass
                else:
                    s.sendall(data)
            else:
                data = s.recv(4096)
                if not data:
                    os._exit(0)
                os.write(1, data)

if __name__ == "__main__":
    main()
