
;; you need to install
;; tabbar
;; yasnippet, ESS, polymode, helm, toc-org, fill-column-indicator
;; paredit
;; other things: ivy+counsel+swiper

(require 'package)
(add-to-list 'package-archives '("melpa" . "https://melpa.org/packages/") t)
;; Comment/uncomment this line to enable MELPA Stable if desired.  See `package-archive-priorities`
;; and `package-pinned-packages`. Most users will not need or want to do this.
;;(add-to-list 'package-archives '("melpa-stable" . "https://stable.melpa.org/packages/") t)
(package-initialize)

(server-start)

;; (when (>= emacs-major-version 24)
;;   (require 'package)
;; ;;  (add-to-list 'package-archives '("melpa" . "http://melpa.milkbox.net/packages/") t)
;; ;;  (add-to-list 'package-archives '("marmalade" . "http://marmalade-repo.org/packages/"))
;; 	;; (add-to-list 'package-archives '("ELPA", "http://tromey.com/elpa/"))
;;   (package-initialize)
;;   )

;; not very useful actually
;; (electric-pair-mode 1)
;; (setq electric-pair-inhibit-predicate
;;       (lambda (c)
;;         (if (char-equal c ?\<) t (electric-pair-default-inhibit c))))

;;;;;;;;;;;;;;;;;;;;
;;;;;;; fonts ;;;;;;
;;;;;;;;;;;;;;;;;;;;
;; I must say the default emacs font is pretty good to the eyes already
;; (add-to-list 'default-frame-alist '(font . "Inconsolata-13"))
;; (add-to-list 'default-frame-alist '(font . "Anonymous Pro-13"))
;; (add-to-list 'default-frame-alist '(font . "Source Code Pro 10"))
;; (set-default-font "Source Code Pro 10")
(add-to-list 'default-frame-alist '(font . "Source Code Pro-13"))
;; (add-to-list 'default-frame-alist '(font . "Consolas-12"))
;; (add-to-list 'default-frame-alist '(font . "Latin Modern Mono-12"))
;; (add-to-list 'default-frame-alist '(font . "Monaco-12"))
;; (add-to-list 'default-frame-alist '(font . "Iosevka-12:light"))

;; auto-mode-alist (append (list '("\\.c$" . c-mode)	
;; 			      '("\\.tex$" . latex-mode)
;; 			      '("\\.S$" . S-mode)
;; 			      '("\\.s$" . S-mode)
;; 			      '("\\.R$" . R-mode)
;; 			      '("\\.r$" . R-mode)
;; 			      '("\\.html$" . html-mode)
;;             '("\\.emacs" . emacs-lisp-mode)
;; 						)
;; 						auto-mode-alist)

;; for Rmarkdown to find pandoc
;; (setq markdown-command "/Applications/RStudio.app/Contents/MacOS/pandoc")
(setq markdown-command "/usr/local/bin/pandoc")

;; open file in existing frame
(setq ns-pop-up-frames nil)

(setenv "PATH" (concat (getenv "PATH") ":/usr/local/texlive/2014/bin/x86_64-darwin/"))
(setq exec-path (append exec-path '("/usr/local/texlive/2014/bin/x86_64-darwin/")))
;; (load-theme 'adwaita t)
;; (load-theme 'ample-zen t)
;; (load-theme 'airline-badwolf t)
;; (load-theme 'tsdh-light t)
(load-theme 'leuven t)
;; (load-theme 'leuven-dark t)
;; (load-theme 'zenburn t)
;; (load-theme 'misterioso t)
;;  (load-theme 'spacemacs-dark t)
;; (load-theme 'deeper-blue t)
;; (load-theme 'dracula t)
;; (setq-default cursor-type 'bar)

(setq-default indent-tabs-mode nil)

(tabbar-mode)

(yas-global-mode 1)
(add-to-list 'yas-snippet-dirs "~/repo/wzlib/emacs/emacs_linux/wanding-config/snippets")
(yas-reload-all)

(defun switch-to-existing-buffer-other-window (part)
  "Switch to buffer with PART in its name."
  (interactive
   (list (read-buffer-to-switch "Switch to buffer in other window: ")))
  (let ((candidates
     (cl-remove
      nil
      (mapcar (lambda (buf)
            (let ((pos (string-match part (buffer-name buf))))
              (when pos
            (cons pos buf))))
          (buffer-list)))))
    (unless candidates
      (user-error "There is no buffers with %S in its name." part))
    (setq candidates (cl-sort candidates #'< :key 'car))
    (switch-to-buffer-other-window (cdr (car candidates)))))

(defun update-r-session ()
  "update the R session without moving focus"
  (interactive)
  (save-excursion
    (progn
      (switch-to-existing-buffer-other-window "*R")
      (end-of-buffer)
      (select-window (previous-window)))))

(global-set-key (kbd "<f10>") `update-r-session)

(setq tab-width 2)
(setq-default tab-width 2)
(setq sh-basic-offset 2
      sh-indentation 2)
(setq-default indent-tabs-mode nil)
(setq org-startup-indented t)
(setq org-startup-folded t)
;; (setq org-startup-folded 'content)
(setq org-src-fontify-natively t)
(setq org-edit-src-content-indentation 0)
(setq org-src-tab-acts-natively t)

(setq org-todo-keywords
      '((sequence "TODO"
                  "|"
                  "DONE"
                  "OBSOLETE"
                  "DEFERRED")))

(if (require 'toc-org nil t)
    (add-hook 'org-mode-hook 'toc-org-mode)
  (warn "toc-org not found"))

(setq mouse-wheel-scroll-amount '(2))
(setq mouse-wheel-progressive-speed nil)
(set-default 'truncate-lines nil)
(global-set-key (kbd "M-q") 'toggle-truncate-lines)

(setq word-wrap nil)
;; command as control
(setq mac-command-modifier 'control)

(add-to-list 'default-frame-alist '(height . 50))
(add-to-list 'default-frame-alist '(width . 100))
(setq-default fill-column 79)
(setq-default c-default-style "k&r"
              c-basic-offset 2)
(setq vc-follow-symlinks t) ;; always follow symlinks when opening file

(show-paren-mode 1)
(require 'mouse)
(xterm-mouse-mode t)
(defun track-mouse (e)) 
(setq mouse-sel-mode t)
(delete-selection-mode 1)

(global-set-key [mouse-4] (lambda () (interactive) (scroll-down 1)))
(global-set-key [mouse-5] (lambda () (interactive) (scroll-up 1)))

(global-set-key (kbd "<f1>") 'ess-display-help-on-object)
;; step by function or region
(global-set-key (kbd "<f12>") 'ess-eval-region-or-function-or-paragraph-and-step)
;; step by line
(global-set-key (kbd "<f11>") 'ess-eval-region-or-line-visibly-and-step)

(defun rename-file-and-buffer (new-name)
  ;; source: http://steve.yegge.googlepages.com/my-dot-emacs-file
  "Renames both current buffer and file it's visiting to NEW-NAME."
  (interactive "sNew name: ")
  (let ((name (buffer-name))
    (filename (buffer-file-name)))
    (if (not filename)
    (message "Buffer '%s' is not visiting a file!" name)
      (if (get-buffer new-name)
      (message "A buffer named '%s' already exists!" new-name)
    (progn
      (rename-file name new-name 1)
      (rename-buffer new-name)
      (set-visited-file-name new-name)
      (set-buffer-modified-p nil))))))

(defun move-buffer-file (dir)
  "Moves both current buffer and file it's visiting to DIR." 
  (interactive "DNew directory: ")
  (let* ((name (buffer-name))
	 (filename (buffer-file-name))
	 (dir
	  (if (string-match dir "\\(?:/\\|\\\\)$")
	      (substring dir 0 -1) dir))
	 (newname (concat dir "/" name)))
    
    (if (not filename)
	(message "Buffer '%s' is not visiting a file!" name)
      (progn 
	(copy-file filename newname 1) 
	(delete-file filename)
 	(set-visited-file-name newname)
 	(set-buffer-modified-p nil)
 	t)))) 

(defun delete-this-buffer-and-file ()
  "Removes file connected to current buffer and kills buffer.
   with dired-x just do 'C-x C-j D yes RET y RET"
  (interactive)
  (let ((filename (buffer-file-name))
        (buffer (current-buffer))
        (name (buffer-name)))
    (if (not (and filename (file-exists-p filename)))
        (error "Buffer '%s' is not visiting a file!" name)
      (when (yes-or-no-p "Are you sure you want to remove this file? ")
        (delete-file filename)
        (kill-buffer buffer)
        (message "File '%s' successfully removed" filename)))))

(defun jao-copy-line ()
  "Copy current line in the kill ring"
  (interactive)
  (kill-ring-save (line-beginning-position)
		  (line-beginning-position 2))
  (message "line copied"))
(global-set-key (kbd "M-l") 'jao-copy-line)

;; M-x print-to-pdf
(defun print-to-pdf ()
  (interactive)
  (ps-spool-buffer-with-faces)
  (switch-to-buffer "*PostScript*")
  (write-file "/tmp/tmp.ps")
  (kill-buffer "tmp.ps")
  (setq cmd (concat "ps2pdf14 /tmp/tmp.ps " (buffer-name) ".pdf"))
  (shell-command cmd)
  (shell-command "rm /tmp/tmp.ps")
  (message (concat "Saved to:  " (buffer-name) ".pdf"))
  )

(defun move-line-upward ()
  (interactive)
  (kill-region (line-beginning-position) (line-beginning-position 2))
  (goto-char (line-beginning-position 0))
  (yank)
  (goto-char (line-beginning-position 0))
  )


(defun move-line-downward ()
  (interactive)
  (kill-region (line-beginning-position) (line-beginning-position 2))
  (goto-char (line-beginning-position 2))
  (yank)
  (goto-char (line-beginning-position 0))
)

(defun insert-before-line ()
  "insert the content yanked (e.g., copied) before the current line"
  (interactive)
  (let ((pos (point))
        (cur-max (point-max)))
    (beginning-of-line)

    ;; I've changed the order of (yank) and (indent-according-to-mode)
    ;; in order to handle the case when yanked line comes with its own indent
    (yank);; (flush-blank-lines (point))(indent-according-to-mode)
    ;; ;; could be as well changed to simple (newline) it's matter of taste
    ;; ;; and of usage
    ;; (newline-and-indent) 
    ;; back to the original position
    ;; (goto-char (+ pos (- (point-max) cur-max)))
    )
  )
(global-set-key (kbd "C-.") 'insert-before-line)

(defun open-line-above ()
  "Open a line above the line the point is at.
Then move to that line and indent according to mode"
  (interactive)
  (move-beginning-of-line 1)
  (newline)
  (previous-line)
  (indent-according-to-mode))
(global-set-key (kbd "C-o") 'open-line-above)
(global-set-key (kbd "C-<up>") 'open-line-above)

(defun open-line-below ()
  "Open a line below the line the point is at.
Then move to that line and indent according to mode"
  (interactive)
  (move-end-of-line 1)
  (newline)
  (indent-according-to-mode))
(global-set-key (kbd "C-,") 'open-line-below)
;; (global-set-key (kbd "C-m") 'open-line-below)

;; (global-set-key (kbd "C-<down>") 'open-line-below)

(defun reload-dot-emacs()
  "Reload .emacs on the fly"
  (interactive)
  (if(bufferp (get-file-buffer ".emacs"))
      (save-buffer(get-buffer ".emacs")))
  (load-file "~/.emacs")
  (message ".emacs reloaded successfully"))

(defun server-shutdown ()
  "Save buffers, Quit, and Shutdown (kill) server"
  (interactive)
  (save-some-buffers)
  (kill-emacs)
  )

(defun my-put-file-name-on-clipboard ()
  "Put the current file name on the clipboard"
  (interactive)
  (let ((filename (if (equal major-mode 'dired-mode)
                      default-directory
                    (buffer-file-name))))
    (when filename
      (let ((x-select-enable-clipboard t)) (kill-new filename))
      ;; (with-temp-buffer
      ;;   (insert filename)
      ;;   (clipboard-kill-region (point-min) (point-max)))
      (message filename))))

(defun comment-or-uncomment-region-or-line ()
    "Comments or uncomments the region or the current line if there's no active region."
    (interactive)
    (let (beg end)
        (if (region-active-p)
            (setq beg (region-beginning) end (region-end))
            (setq beg (line-beginning-position) end (line-end-position)))
        (comment-or-uncomment-region beg end);))
        (next-line)))
(global-set-key (kbd "C-;") 'comment-or-uncomment-region-or-line)

(defun insert-date ()
  "Insert current date yyyy-mm-dd."
  (interactive)
  (when (region-active-p)
    (delete-region (region-beginning) (region-end) )
    )
  (insert (format-time-string "%Y%m%d_"))
  )

(global-set-key "\C-x\C-d" 'insert-date)

(defun timestamp ()
  (interactive)
  (insert (format-time-string "%Y-%m-%d %H:%M")))

;;;;;;;; UNCLEAR FUNCTIONS ;;;;;;;;;;
(defun insert-and-indent-line-above ()
  "[UNCLEAR FUNC] insert the content before the current line and indent the current line"
  (interactive)
  (push-mark)
  (let* 
    ((ipt (progn (back-to-indentation) (point)))
     (bol (progn (move-beginning-of-line 1) (point)))
     (indent (buffer-substring bol ipt)))
    (newline)
    (previous-line)
    (insert indent)))

;; ;;Place all backup copies of files in a common location
;; (defconst use-backup-dir t)   
;; (setq backup-directory-alist (quote ((".*" . "~/emacs-meta/backups/")))
;;       version-control t                ; Use version numbers for backups
;;       kept-new-versions 16             ; Number of newest versions to keep
;;       kept-old-versions 2              ; Number of oldest versions to keep
;;       delete-old-versions t            ; Ask to delete excess backup versions?
;;       backup-by-copying-when-linked t) ; Copy linked files, don't rename.

(global-set-key (kbd "M-i") 'kill-whole-line)
(global-set-key (kbd "<f8>") 'switch-to-buffer)
(global-set-key (kbd "<f9>") 'other-window)
(global-set-key (kbd "<f18>") 'delete-other-windows)
(global-set-key (kbd "<clear>") 'goto-line)
(global-set-key (kbd "M-<up>") 'move-line-upward)
(global-set-key (kbd "M-<down>") 'move-line-downward)
(global-set-key (kbd "<home>") 'beginning-of-line)
(global-set-key (kbd "<end>") 'end-of-line)
(global-set-key (kbd "C-c C-q") 'fill-paragraph)

;; helm
(setq recentf-max-saved-items 10000)
(global-set-key (kbd "M-r") 'helm-recentf)

(defun insert-file-name ()
  "Insert the full path file name into the current buffer."
  (interactive)
  (insert (buffer-file-name (window-buffer (minibuffer-selected-window)))))
(put 'upcase-region 'disabled nil)

;;;;;;;;;;
;; ESS
;;;;;;;;;;
(require 'ess-r-mode)

;; YOU NEED TO UPDATE R VERSION FOR EVERY UPDATE
;; (setq inferior-R-program-name "/Users/zhouw3/.Renv/versions/3.6.1/bin/R")
(setq inferior-R-program-name "/Users/zhouw3/.Renv/versions/4.1.dev/bin/R")


;; Note that there are a lot of out-dated information out
;; there on the internet. ESS has changed a lot in terms of
;; the variables for customization. Look at ess-custom.el
;; and C-h v ess-style-alist for the most updated information.
;;
;; Make sure you use the most updated version of ESS!!
;;
;; To customize, you need to change the default list to OWN
;; all the other style are FIXED!!
;; The default style is RRR. Even changing this thing to
;; DEFAULT won't work. So far only the following works.
(custom-set-variables
 ;; custom-set-variables was added by Custom.
 ;; If you edit it by hand, you could mess it up, so be careful.
 ;; Your init file should contain only one such instance.
 ;; If there is more than one, they won't work right.
 '(custom-safe-themes
   '("bffa9739ce0752a37d9b1eee78fc00ba159748f50dc328af4be661484848e476" "da53c5d117ebada2664048682a345935caf8e54094a58febd5f021462ef20ba2" "614a8fc7db02cb99d9f1acf1297b26f8224cf80bf6c0ec31d30c431503e8b59f" "f2c35f8562f6a1e5b3f4c543d5ff8f24100fae1da29aeb1864bbc17758f52b70" "fa2af0c40576f3bde32290d7f4e7aa865eb6bf7ebe31eb9e37c32aa6f4ae8d10" default))
 '(ess-own-style-list
   '((ess-indent-offset . 4)
     (ess-offset-arguments . prev-line)
     (ess-offset-arguments-newline . prev-line)
     (ess-offset-block . prev-line)
     (ess-offset-continued . straight)
     (ess-align-nested-calls)
     (ess-align-arguments-in-calls)
     (ess-align-continuations-in-calls)
     (ess-align-blocks control-flow)
     (ess-indent-from-lhs)
     (ess-indent-from-chain-start)
     (ess-indent-with-fancy-comments . t)))
 '(package-selected-packages
   '(paredit expand-region doom-themes spacemacs-theme ess laguna-theme zencoding-mode anti-zenburn-theme hc-zenburn-theme labburn-theme zenburn-theme leuven-theme org-drill projectile abs-mode markdown-mode yaml-mode yasnippet-snippets yasnippet toc-org tabbar helm fill-column-indicator)))

(setq-default ess-indent-offset 4)

(setq ess-style 'OWN)
(ess-toggle-underscore nil)

;; (defun rmd-mode ()
;;  "ESS Markdown mode for rmd files"
;;  (interactive)
;;  (require 'poly-R)
;;  (require 'poly-markdown)     
;;  (poly-markdown+r-mode))
(add-to-list 'auto-mode-alist '("\\.Rmd" . poly-markdown+r-mode))
;; the shortcut for R
(defun then_R_operator ()
  "R - %>% operator or 'then' pipe operator"
  (interactive)
  (just-one-space 1)
  (insert "%>%")
  ;; (reindent-then-newline-and-indent)
  (just-one-space 1))
(global-set-key (kbd "C-%") 'then_R_operator)
;; (define-key ess-mode-map (kbd "C-%") 'then_R_operator)
;; (define-key inferior-ess-mode-map (kbd "C-%") 'then_R_operator)

(setq column-number-mode t)

(defun create-tags (dir-name)
     "Create tags file."
     (interactive "DDirectory: ")
     (eshell-command 
      (format "find %s -type f -name \"*.[ch]\" | etags -" dir-name)))

(defun zwd-clean-path ()
  "Clean path to ~ form"
  (interactive)
  (save-restriction
    (narrow-to-region (line-beginning-position) (line-end-position))
    (goto-char (point-min))
    (while (search-forward "/Users/zhouw3" nil t)
      (replace-match "~"))))

(global-set-key (kbd "M-`") 'zwd-clean-path)

(defun xah-select-text-in-quote ()
  "Select text between the nearest left and right delimiters.
Delimiters here includes the following chars: '\"`<>(){}[]“”‘’‹›«»「」『』【】〖〗《》〈〉〔〕（）
This command select between any bracket chars, not the inner text of a bracket. For example, if text is

 (a(b)c▮)

 the selected char is “c”, not “a(b)c”.

URL `http://ergoemacs.org/emacs/modernization_mark-word.html'
Version 2020-03-11"
  (interactive)
  (let (
        ($skipChars "^'\"`<>(){}[]“”‘’‹›«»「」『』【】〖〗《》〈〉〔〕（）〘〙")
        $p1
        )
    (skip-chars-backward $skipChars)
    (setq $p1 (point))
    (skip-chars-forward $skipChars)
    (set-mark $p1)))

(defun xah-select-line ()
  "Select current line. If region is active, extend selection downward by line.
URL `http://ergoemacs.org/emacs/modernization_mark-word.html'
Version 2017-11-01"
  (interactive)
  (if (region-active-p)
      (progn
        (forward-line 1)
        (end-of-line))
    (progn
      (end-of-line)
      (set-mark (line-beginning-position)))))

(defun xah-select-block ()
  "Select the current/next block of text between blank lines.
If region is active, extend selection downward by block.

URL `http://ergoemacs.org/emacs/modernization_mark-word.html'
Version 2019-12-26"
  (interactive)
  (if (region-active-p)
      (re-search-forward "\n[ \t]*\n" nil "move")
    (progn
      (skip-chars-forward " \n\t")
      (when (re-search-backward "\n[ \t]*\n" nil "move")
        (re-search-forward "\n[ \t]*\n"))
      (push-mark (point) t t)
      (re-search-forward "\n[ \t]*\n" nil "move"))))

(defun zwd-path-mark-one ()
  (interactive)
  (if (looking-at "[-_A-Za-z0-9~]")
      (progn
        (while (looking-at "[-_A-Za-z0-9~]") (backward-char))
        (forward-char)
        (set-mark (point))
        (while (looking-at "[-_A-Za-z0-9~]") (forward-char))
        )))

(defun zwd-path-expand-right ()
  (interactive)
  (if (region-active-p)
      (progn
        (if (> (mark) (point)) (exchange-point-and-mark))
        (goto-char (region-end))
        (forward-char)
        (while (looking-at "[-_A-Za-z0-9~]") (forward-char)))))

(defun zwd-path-expand-left ()
  (interactive)
  (if (region-active-p)
      (progn
        (if (< (mark) (point)) (exchange-point-and-mark))
        (goto-char (region-beginning))
        (backward-char)
        (while (looking-back "[-_A-Za-z0-9~]") (backward-char)))))

(defun zwd-path-shrink-left ()
  (interactive)
  (if (region-active-p)
      (progn
        (if (< (mark) (point)) (exchange-point-and-mark))
        (goto-char (region-beginning))
        (while (looking-at "[-_A-Za-z0-9~]") (forward-char))
        (forward-char))))

(defun zwd-path-shrink-right ()
  (interactive)
  (if (region-active-p)
      (progn
        (if (> (mark) (point)) (exchange-point-and-mark))
        (goto-char (region-end))
        (while (looking-back "[-_A-Za-z0-9~]") (backward-char))
        (backward-char))))

(require 'expand-region)
(global-set-key (kbd "M-1") 'zwd-path-expand-left)
(global-set-key (kbd "M-2") 'zwd-path-shrink-left)
(global-set-key (kbd "M-3") 'zwd-path-shrink-right)
(global-set-key (kbd "M-4") 'zwd-path-expand-right)
;; from the biggest (paragraph) to the smallest (one word)
(global-set-key (kbd "M-6") 'er/mark-paragraph)
(global-set-key (kbd "M-7") 'xah-select-line)
(global-set-key (kbd "M-8") 'xah-select-text-in-quote)
(global-set-key (kbd "M-9") 'zwd-path-mark-one)

(custom-set-faces
 ;; custom-set-faces was added by Custom.
 ;; If you edit it by hand, you could mess it up, so be careful.
 ;; Your init file should contain only one such instance.
 ;; If there is more than one, they won't work right.
 )


