from collections import namedtuple

try:
    from psims.controlled_vocabulary.controlled_vocabulary import load_psims
    cv_psims = load_psims()
except Exception:  # pragma: no cover
    cv_psims = None


class Term(namedtuple("Term", ("name", "id", "category", "specialization"))):

    def __eq__(self, other):
        if isinstance(other, str):
            return self.name == other or self.id == other
        else:
            return tuple(self) == tuple(other)

    def __str__(self):
        return self.name

    def __repr__(self):
        text = super(Term, self).__repr__()
        return self.__class__.__name__ + text[4:]

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash(self.name)

    def is_a(self, term):
        """Test whether this entity is exactly ``term`` or a specialization
        of ``term``

        Parameters
        ----------
        term : str or :class:`Term`
            The entity to compare to

        Returns
        -------
        bool
        """
        return term == self.name or term in self.specialization


def _unique_list(items):  # pragma: no cover
    seen = set()
    out = []
    for x in items:
        if x in seen:
            continue
        seen.add(x)
        out.append(x)
    return out


def type_path(term, seed):  # pragma: no cover
    path = []
    i = 0
    steps = []
    try:
        steps.append(term.is_a.comment)
    except AttributeError:
        steps.extend(t.comment for t in term.is_a)
    except KeyError:
        pass
    while i < len(steps):
        step = steps[i]
        i += 1
        path.append(step)
        term = cv_psims[step]
        try:
            steps.append(term.is_a.comment)
        except AttributeError:
            steps.extend(t.comment for t in term.is_a)
        except KeyError:
            continue
    return _unique_list(path)


def render_list(seed, list_name=None, term_cls_name="Term"):  # pragma: no cover
    component_type_list = [seed]
    i = 0
    seen = set()
    if list_name is None:
        list_name = seed.replace(" ", "_") + 's'
    print("%s = [" % (list_name,))
    while i < len(component_type_list):
        component_type = component_type_list[i]
        i += 1
        for term in cv_psims[component_type].children:
            if term.name in seen:
                continue
            seen.add(term.name)
            print("    %s(%r, %r, %r, %r), " % (
                term_cls_name, term.name, term.id, component_type_list[0], type_path(term, seed)))
            if term.children:
                component_type_list.append(term.name)
    print("]")
