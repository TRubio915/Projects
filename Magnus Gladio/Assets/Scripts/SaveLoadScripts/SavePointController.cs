using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SavePointController : MonoBehaviour
{

    public float radiusToCheckForPlayer;
    private bool playerInRange;
    float fade = 0;
    Material mat;
    public LayerMask whatIsPlayer;

    public GameObject savePointText;

    // Start is called before the first frame update
    void Start()
    {
        mat = GetComponent<SpriteRenderer>().material;
        mat.SetFloat("_Fade", fade);
        savePointText.SetActive(false);
    }

    // Update is called once per frame
    void Update()
    {

        playerInRange = Physics2D.OverlapCircle(transform.position, radiusToCheckForPlayer, whatIsPlayer);


        if(playerInRange)
        {
            if (fade <= 0)
            {
                StartCoroutine(EnableFade(1));
                savePointText.SetActive(true);
            }
        }
        else
        {
            if (fade >= 1)
            {
                StartCoroutine(EnableFade(0));
                savePointText.SetActive(false);
            }
        }

    }

    IEnumerator EnableFade(int on)
    {
        if (on == 1)
        {
            while (fade <= 1)
            {
                yield return null;
                mat.SetFloat("_Fade", fade);

                fade += Time.deltaTime * 0.5f;

            }
        }
        else
        {
            while (fade > 0)
            {
                yield return null;

                mat.SetFloat("_Fade", fade);

                fade -= Time.deltaTime * 0.5f;

            }
        }
    }


    private void OnDrawGizmos()
    {
        Gizmos.color = Color.black;
        Gizmos.DrawWireSphere(transform.position, radiusToCheckForPlayer);
    }
}
