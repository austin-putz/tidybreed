-- Lua filter for tidybreed HTML docs
-- 1. Wrap tables in a scroll div
-- 2. Turn Y / M / N cells into coloured badges

function Table(el)
  return pandoc.Div(el, pandoc.Attr("", {"table-wrap"}))
end

function Str(el)
  local s = el.text
  if s == "Y" then
    return pandoc.RawInline("html", '<span class="badge-Y">Y</span>')
  elseif s == "M" then
    return pandoc.RawInline("html", '<span class="badge-M">M</span>')
  elseif s == "N" then
    return pandoc.RawInline("html", '<span class="badge-N">N</span>')
  end
  return el
end
