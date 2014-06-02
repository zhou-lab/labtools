(require 'w3m-search)
(eval-after-load "w3m-search"
  '(progn
     (add-to-list 'w3m-search-engine-alist
		  '("dictcn"
		    ;; "http://www.google.com.hk/dictionary?langpair=en|zh-CN&hl=zh-CN&q=%s"
		    ;; "http://dict.cn/%s"
		    "http://hk.dictionary.yahoo.com/dictionary?p=%s"
		    nil))
     (add-to-list 'w3m-search-engine-alist
		  '("freedict"
		    ;; "http://www.google.com.hk/dictionary?langpair=en|zh-CN&hl=zh-CN&q=%s"
		    "http://www.thefreedictionary.com/%s"
		    nil))
     (add-to-list 'w3m-uri-replace-alist
		  '("\\'dc:" w3m-search-uri-replace "dictcn"))))

;; w3m-search
(defun search-chinese (query)
  "search dict.cn for chinese translation of an english word"
  (interactive "sEnter the english word: ")
  (w3m-search "dictcn" query))
(global-set-key (kbd "<f6>") 'search-chinese)

(defun search-english (query)
  "search thefreedictionary.com for english translation of an english word"
  (interactive "sEnter the english word: ")
  (w3m-search "freedict" query))
(global-set-key (kbd "<f7>") 'search-english)

(defun search-google (query)
  "search google"
  (interactive "sEnter the search string: ")
  (w3m-search "google" query))


(require 'w3m)
(defun w3m-add-keys ()
  (define-key w3m-mode-map "\M-r" 'anything-recentf)
  (define-key w3m-mode-map "\M-b" 'anything-buffers+)
  (define-key w3m-mode-map "\M-f" 'anything-find-files))
(add-hook 'w3m-mode-hook 'w3m-add-keys)
(setq w3m-use-cookies t)
;; (require 'w3m-search)
;; (add-to-list 'w3m-search-engine-alist
;; 	     '("emacs-wiki" "http://www.emacswiki.org/cgi-bin/wiki.pl?search=%s"))
