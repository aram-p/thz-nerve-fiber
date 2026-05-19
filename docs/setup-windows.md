# Windows setup — handoff to Claude Code

This file is the one-stop walkthrough for getting the project running on your Windows laptop. The first half is ~10 minutes of manual setup you do once. The second half is a prompt you hand to Claude Code, after which it works through the GitHub issues on its own.

---

## Part 1 — Bootstrap (you, ~10 minutes)

Open **PowerShell as Administrator** (right-click Start → "Terminal (Admin)" or "Windows PowerShell (Admin)"). Run these in order.

### 1. Install everything in one go

```powershell
winget install tailscale.tailscale
winget install Git.Git
winget install GitHub.cli
winget install astral-sh.uv
winget install Anthropic.Claude
Add-WindowsCapability -Online -Name OpenSSH.Server~~~~0.0.1.0
Start-Service sshd
Set-Service -Name sshd -StartupType Automatic
```

You may need to press `Y` a few times to accept package terms. After this finishes, **close PowerShell and open a fresh one** (still as Administrator) so the new tools land on PATH.

### 2. Log in to Tailscale

```powershell
tailscale up
```

A browser opens. Sign in with whatever account you used on the Mac (Google / GitHub / Microsoft). Once it says "Success," close the browser tab.

### 3. Log in to GitHub

```powershell
gh auth login
```

Pick: **GitHub.com** → **HTTPS** → **Y** (authenticate Git with GitHub credentials) → **Login with a web browser** → copy the one-time code shown → press Enter → paste the code in the browser → done.

### 4. Clone the repo

```powershell
cd $env:USERPROFILE
gh repo clone aram-p/thz-nerve-fiber
cd thz-nerve-fiber
```

### 5. Start Claude Code

```powershell
claude
```

A Claude Code session opens inside the repo.

---

## Part 2 — Paste this prompt into Claude

Copy everything inside the box below and paste it as your first message to Claude. (Triple-click any line to select that line; or use the GitHub "Copy" button if you're reading this on github.com.)

```
You are Claude Code running on a Windows laptop, inside the cloned `thz-nerve-fiber` repository. Your job is to drive this thesis simulation project autonomously.

== Orientation (do this first) ==

1. Read `CLAUDE.md` and every file under `docs/adr/`. These are user-approved decisions; do not modify them.
2. List open GitHub issues: `gh issue list --state open --limit 30`. Read issues #2 through #14 in full: `gh issue view <number>`.
3. Verify tooling: `uv --version`, `git --version`, `gh --version`. They should all print versions.
4. Locate `comsolbatch.exe`: `Get-ChildItem "C:\Program Files\COMSOL" -Recurse -Filter comsolbatch.exe -ErrorAction SilentlyContinue | Select-Object FullName`. Note the path; you'll need it.

== Workflow ==

Work through the issues in numeric order, skipping any labeled `ready-for-human`. For each issue:

1. `gh issue view <N>` — read the body fully.
2. Plan the work briefly (one paragraph) — confirm the plan to the user only if it deviates from the issue body.
3. Implement.
4. Run the "Done when" check from the issue body. If it fails, debug once. If it still fails, stop and ask.
5. Commit with message format: `Phase X.Y: <short summary>` and a body that says `Closes #N`. One commit per issue.
6. `gh issue close <N> --comment "<brief summary of what you did + any caveats>"`.
7. `git push`.
8. Move to the next open issue.

== Stop and ask the user when ==

- `comsolbatch.exe` is not found on the system
- Any test or "Done when" check fails twice
- An issue is labeled `ready-for-human` (skip it, note the skip, continue)
- You're about to do anything destructive (`rm -rf`, force-push, drop a database, uninstall something)
- You hit an actual design choice that wasn't pre-decided in CLAUDE.md or the ADRs

== Do not ==

- Modify ADRs unilaterally — propose changes, don't enact them
- Skip the validation issue (#11) and act like Phase 2 is done — it is gated on the user's thesis advisor
- Open the COMSOL GUI; everything runs via `mph` and `comsolbatch`

== Begin ==

Start with the orientation steps. When you've finished those and confirmed all four tools work, post a one-line status ("Orientation done, beginning issue #2") and proceed.
```

---

## Part 3 — Finishing touches (do later, when convenient)

These aren't blocking — Claude can start working without them. But you'll want them before you leave the house with your Mac.

### Power settings (so the laptop stays awake)

- Settings → System → Power & battery → **Screen and sleep** → "When plugged in, put my device to sleep after" → **Never**
- Control Panel → Power Options → "Choose what closing the lid does" → "When I close the lid, plugged in" → **Do nothing**
- Device Manager → Network adapters → your Wi-Fi card → Properties → **Power Management** tab → uncheck "Allow the computer to turn off this device to save power"

### Tailscale rename

Open https://login.tailscale.com/admin/machines on any device. Find your Windows machine in the list, `...` menu → "Edit machine name" → rename to `thz-tower`. Confirm **MagicDNS** is on.

### SSH key auth (passwordless from Mac)

From your **Mac**:

```
ssh-copy-id <your-windows-username>@thz-tower
```

If that fails (Windows OpenSSH is fussy), manually paste the contents of `~/.ssh/id_*.pub` (from Mac) into `C:\Users\<user>\.ssh\authorized_keys` on Windows.

Test: `ssh <user>@thz-tower` should connect without a password prompt.

---

## What success looks like

When Part 1 is done, you have a Windows laptop with Claude Code running inside the repo. Paste the Part 2 prompt and walk away. When you come back:

- Issues #2–#10 (except #11) should be closed
- The `main` branch has ~9 new commits
- The Python pipeline is ready to run baseline sweeps
- Claude has stopped at #11 awaiting your validation call (probably with your advisor)

At that point you read the diff, run the baseline yourself, and decide whether to greenlight Phase 3 or rework.
